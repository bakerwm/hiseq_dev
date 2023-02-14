#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Demultiplex fastq by barcode or index; (single)
input:
  - fq1 : read1 of paired-end
  - fq2 : read2 of paired-end, or None
  - out_dir : directory saving files
  - barcode.csv : format "name,seq"
  - mismatch : number of mismatches allowed
  - barcode_n_left: int
  - barcode_n_right: int
  - parallel_jobs : int
output:
  - out_dir/report.txt
  - out_dir/read_count.json
  - out_dir/...
Function:
1. Demultiplex index, single: '2:N:0:ATCACGAT+AGATCTCG'; i7->i5
    - 1:N:0:ATCACGAT+AGATCTCG (i7+i5)
    - 1:N:0:ATCACGAT (i7)
## others
1. defq, ultra fasta multi-threaded fq demultiplexing
https://github.com/OpenGene/defq
2. deML, maxlikelihood demultiplexing
https://github.com/grenaud/deML
## test
step1: split by i7-index: index1_index2.r1.fq
step2: split by inline barcode:
## search index
index1:index2:barcode
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from itertools import combinations
from contextlib import ExitStack
import Levenshtein as lev # distance
from hiseq.utils.utils import log, update_obj, Config, get_date, str_distance
from hiseq.utils.file import check_dir, file_abspath, file_exists
from hiseq.utils.seq import Fastx, check_fx_paired, check_fx_args, readfq


class DemxR1(object):
    """
    Demultiplex barcode, at the beginning of read1/2
    Example:
    >>> args = {
        'fq1': '1m.r1.fq.gz',
        'fq2': '1m.r2.fq.gz',
        'out_dir': 'aaa',
        'index_table': 'index.csv',
        'mismatch': 0,
        'demo': False,
        'gzipped': True,
    }
    >>> Demx1(**args).run()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'fq1': None,
            'fq2': None,
            'index_type': 'i7', # i7, barcode
            'out_dir': None,
            'barcode_in_read': 1,
            'barcode_n_left': 0,
            'barcode_n_right': 1,
            'index_table': None, #name,index
            'mismatch': 0,
            'threads': 1,
            'overwrite': False,
            'demo': False,
            'gzipped': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'demx_r1'
        # fastq file
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq invalid: {}, {}'.format(self.fq1, self.fq2))
        # out_dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        # if self.demo:
        #     self.out_dir = os.path.join(self.out_dir, 'demo')
        if not self.mismatch in range(4):
            raise ValueError('illegal mismatch: [{}], expect [0,1,2,3]'.format(
                self.mismatch))
        # choose demx function
        self.guess_name = self.guess_index_name if self.index_type == 'i7' else self.guess_barcode_name
        # index table
        self.idx = self.load_index(self.index_table, min_cols=2)
        # files
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.index_table = file_abspath(self.index_table)
        self.out_dir = file_abspath(self.out_dir)
        self.fn_json = os.path.join(self.out_dir, 'read_count.json')
        self.report_txt = os.path.join(self.out_dir, 'report.txt')
        check_dir(self.out_dir)


    def load_index(self, x, min_cols=2, parse_read_count=False):
        """
        index format: name, seq
        """
        d = {}
        try:
            with open(x) as r:
                for l in r:
                    if l.startswith('#') or len(l.strip()) == 0:
                        continue
                    s = re.split('[,\s\t]', l.strip())
                    if len(s) < min_cols:
                        raise ValueError(
                            'at least {} cols required: {}'.format(min_cols, x))
                    # if not len(s) == 2:
                    #     raise ValueError('illegal index_table: {}'.format(l))
                    name, i7 = s[:2]
                    if not re.match('^[ACGTN$]+$', i7, flags=re.IGNORECASE):
                        continue
                    if i7 in d:
                        raise ValueError('index not unique, {}'.format(
                            self.index_table))
                    # parse read count
                    p2 = re.compile('^([0-9\.]+)M?$', flags=re.IGNORECASE)
                    g2 = p2.match(s[-1])
                    if parse_read_count:
                        val = eval(g2.group(1)) if g2 else 0
                    else:
                        val = name
                    # fix truseq
                    if self.index_type == 'i7':
                        i7 = self.fix_truseq(i7)
                    d.update({i7:val}) # idx:name
        except Exception as exc:
            print(exc)
        # check i7
        if len(d) == 0:
            raise ValueError('no indexes: {}'.format(self.index_table))
        # check index table
        if not self.is_valid_index(list(d.keys())):
            raise ValueError('index not compatible with mismatch={}'.format(
                self.mismatch))
        # load barcode width
        self.barcode_width = len(list(d.keys())[0]) #
        return d


    def fix_truseq(self, x):
        """
        fix index, mixed with TruSeq (6) and Nextera (8)
        solution: 
        TruSeq: 6 to 8 bases, add "AT" to the tail
        ATCACG -> ATCACGAT
        """
        if isinstance(x, str):
            if len(x) == 6:
                x += 'AT'
        elif isinstance(x, list):
            x = [self.fix_truseq(i) for i in x]
        return x


    def is_valid_index(self, x):
        """
        Parameters:
        -----------
        x : list
            list of index sequences.
        valid index table:
        1. mismatch between any two index
        2. same length ? # no
        """
        # combination for list
        if isinstance(x, list):
            c1 = [str_distance(a, b) for a, b in list(combinations(x, 2))]
            c2 = [i > self.mismatch for i in c1] #
            c3 = len(set(map(len, x))) == 1 # same length
            out = all(c2) and c3
        else:
            log.error('expect list, got {}'.format(type(x).__name__))
            out = False
        return out


    def subset_fq(self, smp='demo', sample_size=1000000):
        """
        Create a subsample of input fastq files, default: 1M reads
        Run the whole process for demonstration
        update: fq1, fq2, out_dir
        """
        log.info('Running demo with subsample: {} reads'.format(sample_size))
        # update args
        self.data_dir = os.path.join(self.out_dir, 'data')
        check_dir(self.data_dir)
        # subset
        demo_fq1 = os.path.join(self.data_dir, os.path.basename(self.fq1))
        Fastx(self.fq1).sample(demo_fq1, sample_size)
        self.fq1 = demo_fq1
        if self.fq2:
            demo_fq2 = os.path.join(self.data_dir, os.path.basename(self.fq2))
            Fastx(self.fq2).sample(demo_fq2, sample_size)
            self.fq2 = demo_fq2


    def guess_index_name(self, x):
        """
        Parameters:
        -----------
        x : list (name, seq, qual, comment)
            comment of fastq file, example: 
            - 1:N:0:ATCACGAT+AGATCTCG (i7+i5)
            - 1:N:0:ATCACGAT (i7)
        Guess the name from index_table
        return "undemx" if not matched
        """
        m = x[-1] # comment
        # extract index seq
        s = re.split('[:+]', m)
        query = s[3] # i7, ignore i5
        # i5 = s[4] if len(s) > 4 else None
        # query = i7 if self.index_type == 'i7' else i5
        # search index name
        hit = None
        for i,j in self.idx.items():
            mm = str_distance(query, i)
            if mm >= 0 and mm <= self.mismatch:
                hit = j
                # To speed up, exit the loop when meeting the first hit.
                # if mm>0, this strategy might skip the best hit, when
                # mismatches hit at first.
                break
        # return the name
        return hit if hit else 'undemx'


    def guess_barcode_name(self, x):
        """
        Parameters:
        -----------
        x : list (name, seq, qual, comment)
            comment of fastq file, example: seq
        Guess the name from index_table
        return "undemx" if not matched
        """
        s = x[1] # seq
        query = s[self.barcode_n_left:(self.barcode_n_left+self.barcode_width)]
        # search index name
        hit = None
        for i,j in self.idx.items():
            mm = str_distance(query, i)
            if mm >= 0 and mm <= self.mismatch:
                hit = j
                # To speed up, exit the loop when meeting the first hit.
                # if mm>0, this strategy might skip the best hit, when
                # mismatches hit at first.
                break
        # return the name
        return hit if hit else 'undemx'


    def demx_se(self, fh):
        # prepare output files, named by index
        fq_names = list(self.idx.values())
        fq_names.append('undemx')
        sx = '.fq.gz' if self.gzipped else '.fq'
        fq_files = [os.path.join(self.out_dir, i + sx) for i in fq_names]
        fn = {}
        with ExitStack() as stack:
            # writers
            fws = [stack.enter_context(xopen(f, 'wt')) for f in fq_files]
            n = 0 # counter
            for r1 in readfq(fh):
                r1 = list(r1) # convert tuple to list
                n += 1
                if n%1e6 == 0:
                    log.info('Processed reads: {}'.format(n))
                if r1[-1] is None:
                    raise ValueError('comment field missing: {}'.format(self.fq1))
                r1[0] = r1[0] + ' ' + r1[-1] # update name
                fq = '\n'.join(['@' + r1[0], r1[1], '+', r1[2]])
                # hit = self.guess_index_name(r1) # seq
                hit = self.guess_name(r1) # seq
                fw = fws[fq_names.index(hit)] # guess writer
                fw.write(fq+'\n')
                fn[hit] = fn.get(hit, 0) + 1
        # save dict to json
        Config().dump(fn, self.fn_json)
        return fq_files


    def demx_pe(self, fh1, fh2):
        # prepare output files, named by index
        fq_names = list(self.idx.values())
        fq_names.append('undemx')
        sx = '.fq.gz' if self.gzipped else '.fq'
        fq1_files = [os.path.join(self.out_dir, i + '_1' + sx) for i in fq_names]
        fq2_files = [os.path.join(self.out_dir, i + '_2' + sx) for i in fq_names]
        fn = {}
        with ExitStack() as stack:
            fws1 = [stack.enter_context(xopen(f, 'wt')) for f in fq1_files]
            fws2 = [stack.enter_context(xopen(f, 'wt')) for f in fq2_files]
            n = 0 # counter
            for r1, r2 in zip(readfq(fh1), readfq(fh2)):
                r1, r2 = [list(r1), list(r2)]
                n += 1
                if n%1e6 == 0:
                    log.info('Processed reads: {}'.format(n))
                if r1[-1] is None:
                    raise ValueError('comment field missing: {}'.format(self.fq1))
                if not r1[0] == r2[0]:
                    raise ValueError('file not paired, check: {}'.format(r1[0]))
                # update name
                r1[0] = r1[0] + ' ' + r1[-1] # comment
                r2[0] = r2[0] + ' ' + r2[-1] # comment
                fq1 = '\n'.join(['@'+r1[0], r1[1], '+', r1[2]])
                fq2 = '\n'.join(['@'+r2[0], r2[1], '+', r2[2]])
                rx = r2 if self.barcode_in_read == 2 else r1 #
                # hit = self.guess_index_name(rx) # seq
                hit = self.guess_name(rx) # seq
                fw1 = fws1[fq_names.index(hit)]
                fw2 = fws2[fq_names.index(hit)]
                fw1.write(fq1 + '\n')
                fw2.write(fq2 + '\n')
                fn[hit] = fn.get(hit, 0) + 1
        # save dict to json
        Config().dump(fn, self.fn_json) # read count
        return [fq1_files, fq2_files]


    def show_msg(self):
        msg = '\n'.join([
            '='*80,
            '{:>20} : {}'.format('Program', 'Demx_index'),
            '{:>20} : {}'.format('Date', get_date()),
            '{:>20} : {}'.format('fq1', self.fq1),
            '{:>20} : {}'.format('fq2', self.fq2 if self.fq2 else 'None'),
            '{:>20} : {}'.format('PE_mode', 'yes' if self.is_paired else 'no'),
            '{:>20} : {}'.format('index_table', self.index_table),
            '{:>20} : {}'.format('out_dir', self.out_dir),
            '{:>20} : {}'.format('report', self.report_txt),
            '{:>20} : {}'.format('index', self.index_type),
            '{:>20} : {}'.format('barcode_in_read', self.barcode_in_read),
            '{:>20} : {}'.format('barcode_width', self.barcode_width),
            '{:>20} : {}'.format('barcode_n_left', self.barcode_n_left),
            '{:>20} : {}'.format('barcode_n_right', self.barcode_n_right),
            '{:>20} : {}'.format('mismatch', self.mismatch),
            '{:>20} : {}'.format('overwrite', 'yes' if self.overwrite else 'no'),
            '{:>20} : {}'.format('demo', 'yes' if self.demo else 'no'),
            '='*80,
        ])
        print(msg)


    def report(self):
        # load expect read count (million)
        df1 = self.load_index(self.index_table, parse_read_count=True)
        total_exp = sum(df1.values())
        # load real read count
        df2 = Config().load(self.fn_json) if file_exists(self.fn_json) else None
        if df2 is None:
            df2 = {}
        total = sum(df2.values())
        if total < 1:
            total = 1000000 # default: 1M
        # check output pct
        output_pct = total/(total_exp*1e4) if total_exp > 0 else 100.0
        i = 0
        f_stat = []
        idx = self.idx.copy() # local
        idx.update({'null':'undemx'})
        for k,v in idx.items():
            i += 1
            n = df2.get(v, 0) # count
            n_exp = df1.get(k, 0) # expect read count (million)
            n_pct =n/(n_exp*1e4) if n_exp > 0 else 100.0
            # elements:
            s = ' '.join([
                '{:>5d}'.format(i),
                '{:<50}'.format(v),
                '{:>12,}'.format(n),
                '{:>8.1f}'.format(n/1e6),
                '{:>8}'.format(n_exp),
                '{:>8.1f}%'.format(n_pct),
                # '{:8.1f}'.format(v/1e6*scale),
            ])
            f_stat.append(s)
        # output
        msg = '\n'.join([
            '='*80,
            '{:>20} : {}'.format('Program', 'Demx (report)'),
            '{:>20} : {}'.format('Date', get_date()),
            '{:>20} : {:6.1f} M'.format('Expect reads', total_exp),
            '{:>20} : {:6.1f} M ({:>11,}) {:>8.1f}%'.format(
                'Total reads', total/1e6, total, output_pct
            ),
            '{:>5} {:<50s} {:>12} {:>8} {:>8} {:>8}'.format(
                'order', 'filename', 'count', 'million', 'expect', 'percent'),
            '\n'.join(f_stat),
            '='*80,
        ])
        # save to file
        with open(self.report_txt, 'wt') as w:
            w.write(msg+'\n')
        print(msg)


    def run(self):
        if file_exists(self.fq2):
            self.is_paired = check_fx_paired(self.fq1, self.fq2)
        else:
            self.is_paired = False
        self.show_msg()
        if file_exists(self.fn_json):
            log.info('demx finished, see file: {}'.format(self.fn_json))
        else:
            if self.demo:
                self.subset_fq()
            if self.is_paired:
                with xopen(self.fq1) as r1, xopen(self.fq2) as r2:
                    self.demx_pe(r1, r2)
            else:
                with xopen(self.fq1) as r1:
                    self.demx_se(r1)
        self.report()


def get_args():
    example = '\n'.join([
        'De-multiplex single index/barcode',
        '1. only for i7 demx',
        '$ python demx1.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t i7',
        '2. only for barcode demx',
        '$ python demx1.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t barcode -x 2 -l 3 -r 2',
    ])
    parser = argparse.ArgumentParser(
        prog='demx1',
        description='De-multiplex single index',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
        help='read1 in fastq format, gzipped')
    parser.add_argument('-2', '--fq2',
        help='read2 in fastq format, gzipped, (optional)')
    parser.add_argument('-s', '--index-table', dest='index_table',
        required=True,
        help='index table in csv format, [filename,barcode]')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    parser.add_argument('-t', '--index-type', dest='index_type', default='i7',
        choices=['i7', 'barcode'],
        help='type of index, [i7, barcode], default: [i7]')
    parser.add_argument('-x', '--barcode-in-read', dest='barcode_in_read',
        choices=[1, 2], default=1, type=int,
        help='barcode in read1/2, default: [1]')
    parser.add_argument('-l', '--barcode-n-left', type=int,
        dest='barcode_n_left', default=0,
        help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int,
        dest='barcode_n_right', default=0,
        help='bases locate on the right of barcode')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-j', '--parallel-jobs', type=int, dest='parallel_jobs',
        default=1, help='number of jobs run in parallel, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    parser.add_argument('-z', '--gzipped', action='store_true',
        help='gzip compress output files')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demonstration, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    DemxR1(**args).run()


if __name__ == '__main__':
    main()

#