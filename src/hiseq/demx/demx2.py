#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1. rename fastq files (fq file were named by: TruSeq_Index1-48; Next_Ad2.1-24)
2. demx bc files
## others
1. defq, ultra fasta multi-threaded fq demultiplexing
https://github.com/OpenGene/defq
2. deML, maxlikelihood demultiplexing
https://github.com/grenaud/deML
"""


import os
# import sys
import re
import pathlib
import argparse
from xopen import xopen
from collections import Counter
from contextlib import ExitStack
from multiprocessing import Pool
from hiseq.demx.demx_r1 import DemxR1
from hiseq.demx.sample_sheet import HiSeqIndex
from hiseq.utils.utils import log, update_obj, Config, get_date
from hiseq.utils.seq import Fastx, list_fx, list_fx2
from hiseq.utils.file import (
    check_dir, file_abspath, file_prefix, file_exists, symlink_file
)


class Demx2(object):
    """
    Demultiplex fastq by i7+barcode (optional)
    Arguments
    ---------
    sample_sheet : str
        The table of sample list, csv file
        ['sample_name', 'i7_seq', 'i5_seq', 'bc_seq', 'reads']
    data_dir, str
        The path to the fastq files
    out_dir, str
        The path to the dir, final output
    Description
    Example:
    >>> args = {
        'fq1': '1m.r1.fq.gz',
        'fq2': '1m.r2.fq.gz',
        'out_dir': 'aaa',
        'sample_sheet': 'index.csv',
        'mismatch': 0,
        'demo': False,
        'gzipped': True,
    }
    >>> Demx2x(**args).run()
    -----------
    The fastq files in data_dir are named by the i7_index_name; in case, some
    of the i7_index file contains multiple sub_files, distinguished by
    in-line barcode
    This function is designed to do:
    1. rename i7_only files, (retrieve sample_name from table, by i7_index_name)
    2. demultiplex the i7 files, contains barcode
    3. organize the report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'sample_sheet': None,
            'data_dir': None,
            'out_dir': None,
            'mismatch': 0,
            'barcode_in_read': 2,
            'barcode_n_left': 0,
            'barcode_n_right': 1,
            'mismatch': 0,
            'overwrite': False,
            'demo': False,
            'gzipped': True,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'demx_r2'
        # out_dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        if not self.mismatch in range(4):
            raise ValueError('illegal mismatch: [{}], expect [0,1,2,3]'.format(
                self.mismatch))
        # check fastq files
        self.raw_fq_list = list_fx(self.data_dir, recursive=True)
        if len(self.raw_fq_list) < 2:
            raise ValueError('no fastq files: {}'.format(self.data_dir))
        # index table
        self.init_files()
        self.idx = self.load_index(self.index_table, min_cols=4)
        self.demx_type, self.index_list = self.split_i7_bc()
        if self.demx_type > 1:
            self.i7_with_bc = file_prefix(self.index_list[1])
            self.i7_with_bc = [i.replace('_table', '') for i in self.i7_with_bc]
        else:
            self.i7_with_bc = [] # empty
        # save config
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_files(self):
        self.sample_sheet = file_abspath(self.sample_sheet)
        self.out_dir = file_abspath(self.out_dir)
        self.config_dir = os.path.join(self.out_dir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        self.fn_json = os.path.join(self.out_dir, 'read_count.json')
        self.i7_fn_json = os.path.join(self.out_dir, 'i7_index', 'read_count.json')
        # self.bc_fn_json = os.path.join(self.out_dir, 'bc_index', 'read_count.json')
        self.report_txt = os.path.join(self.out_dir, 'report.txt')
        self.i7_table = os.path.join(self.out_dir, 'index_table.csv')
        self.bc_table = os.path.join(self.out_dir, 'barcode_table.csv')
        check_dir(self.config_dir)


    def load_index(self, x, min_cols=2, parse_read_count=False):
        """
        index table: name, i7_seq, i5_seq, barcode_seq, ...
        ignore: i5 # !!!
        read_count (Million), at last column
        """
        d = {} # i7+bc
        try:
            with open(x) as r:
                for l in r:
                    if l.startswith('#') or len(l.strip()) == 0:
                        continue
                    s = re.split('[,\s\t]', l.strip()) #
                    if len(s) < min_cols:
                        raise ValueError(
                            'at least {} cols required: {}'.format(min_cols, x))
                    # match format
                    p1 = re.compile('^null$|^[ACGTN]+$', flags=re.IGNORECASE)
                    if len(s) >= 4:
                        name, i7, i5, bc = s[:4]
                        idx = '{}:{}'.format(i7, bc)
                        idx_list = s[2:4]
                    elif len(s) >= 2:
                        name, idx = s[:2]
                        idx_list = [s[1]]
                    else:
                        raise ValueError(
                            'at least {} cols required: {}'.format(2, x))
                    # if len(s) <= 2:
                    #     name, idx = s[:2]
                    #     idx_list = [s[1]]
                    # elif len(s) <= 4:
                    #     name, i7, i5, bc = s[:4]
                    #     idx = '{}:{}'.format(i7, bc)
                    #     idx_list = s[2:4]
                    px = [p1.match(i) is not None for i in idx_list] #
                    if not all(px):
                        continue
                    # parse read count
                    p2 = re.compile('^([0-9\.]+)M?$', flags=re.IGNORECASE)
                    g2 = p2.match(s[-1])
                    if parse_read_count:
                        val = eval(g2.group(1)) if g2 else 0
                    else:
                        val = name
                    d.update({idx:val}) # idx:name
        except Exception as exc:
            log.error(exc)
        # check i7
        if len(d) == 0:
            raise ValueError('no indexes: {}'.format(self.index_table))
        return d


    def split_i7_bc(self):
        """
        split index table: i7->bc->name
        1. into i7
        2. into bc
        3. into i7 + bc
        """
        # sub_dir: index table, named by i7
        p = re.compile('^null$|^[ACGTN]+$', flags=re.IGNORECASE)
        x = [i.split(':') for i in list(self.idx.keys())] # (i7, bc)
        i7 = Counter([i[0] for i in x if p.match(i[0])]) # freq
        bc = Counter([i[1] for i in x if p.match(i[1])]) # freq
        if len(i7) == len(self.idx) or len(bc) == len(self.idx):
            tp = 1 # single
            idx_list = []
            for k,v in self.idx.items():
                k1,k2 = k.split(':') # i7,bc
                idx_list.append(','.join([v, k1]))
            out = self.i7_table if len(i7) == len(self.idx) else self.bc_table
            with open(out, 'wt') as w:
                w.write('\n'.join(idx_list)+'\n')
        else:
            # 3. i7 + bc
            tp = 2 # multi
            ## i7 to N-bc
            i7_dir = os.path.join(self.out_dir, 'i7_index')
            check_dir(i7_dir)
            i7_table = os.path.join(i7_dir, 'i7_table.csv')
            out = [i7_table]
            i7_list = {} # all unique i7
            with open(i7_table, 'wt') as w:
                for k,v in self.idx.items():
                    k1,k2 = k.split(':') # i7,bc
                    if k1 in i7_list:
                        continue
                    if i7.get(k1, 1) > 1:
                        i7v = k1
                        i7_list.update({k1:v})
                    else:
                        i7v = v
                    w.write(','.join([i7v, k1])+'\n')
            ## N-bc
            bc_dir = os.path.join(self.out_dir, 'bc_index')
            check_dir(bc_dir) # all
            bc_files = [
                os.path.join(bc_dir, i+'_table.csv') for i in i7_list
            ]
            out.append(bc_files)
            with ExitStack() as stack:
                fws = [stack.enter_context(xopen(f, 'wt')) for f in bc_files]
                for k,v in self.idx.items():
                    k1,k2 = k.split(':') # i7,bc
                    if i7.get(k1, 1) > 1:
                        bc_f = os.path.join(bc_dir, k1+'_table.csv')
                        fw = fws[bc_files.index(bc_f)]
                        fw.write(','.join([v, k2])+'\n')
        return [tp, out]


    def to_i7_index(self, x):
        """
        index_name: 
        TruSeq_Index13
        Next_Ad2.1
        index_seq:
        Index().()
        """
        # extract i7 name
        p = re.compile('True?Seq[\._-]Index\d{1,2}|Next_Ad2[\._-]\d{1,2}|D7\d+', re.IGNORECASE)
        g = p.search(x)
        if g:
            i7 = g.group() # 
            i7 = re.sub('Ad2.', 'Ad2.', i7, flags=re.IGNORECASE) # fix Ad2[.-_] to Ad2.
            i7 = re.sub('TruSeq.', 'TruSeq_', i7, flags=re.IGNORECASE) # fix TruSeq- to TruSeq_
            # Convert i7 seq
            out = [i7, HiSeqIndex(i7).index]
        else:
            out = None
        return out


    def extract_fq_suffix(self, x):
        """
        Extract the suffix of fastq file:
        - *_R1.fq.gz -> _1.fq.gz (PE)
        - *_2.fq.gz  -> _2.fq.gz (PE)
        - *_r1.fastq -> _1.fq (PE)
        """
        # xname = os.path.basename(x)
        # is_r1 = re.search('_(R)?1.f(ast)?q+.gz', xname, re.IGNORECASE)
        p = re.compile('(\.|_)(R?[12]).(f(ast)?q+)(.gz)?', re.IGNORECASE)
        g = p.search(x) # 'Next_Ad2-1_R1.fq.gz'
        if g:
            s = g.groups() # ('_', 'R1', 'fq', None, '.gz')
            r12 = re.sub('r', '', s[1], flags=re.IGNORECASE)
            out = '_{}.{}'.format(r12, s[2])
            if s[4]:
                out += s[4]
        else:
            out = None
        return out


    def rename_by_i7(self, x):
        """
        fastq file name format:
        - Next_Ad2.1_1.fq.gz
        - Next_Ad2-1_1.fq.gz
        - Next_Ad2_1_1.fq.gz
        - YY130-G1-Next_Ad2_1.1.fq.gz
        - ...
        """
        # 1. single i7 mode
        i7_df = self.load_index(x, min_cols=2) # i7 table
        i7_dir = os.path.join(self.out_dir, 'i7_index')
        i7_fn = os.path.join(i7_dir, 'read_count.json')        
        i7_fn_df = Config().load(i7_fn) if file_exists(i7_fn) else None
        if i7_fn_df is None:
            i7_fn_df = {}
        check_dir(i7_dir)
        # rename files
        for fq in self.raw_fq_list:
            i7 = self.to_i7_index(fq)
            if i7 is None:
                # log.error('unknown fq: {}'.format(fq))
                continue
            i7_id, i7_seq = i7
            name = i7_df.get(i7_seq, None)
            suffix = self.extract_fq_suffix(fq)
            if name is None:
                # log.error('unknown i7: {} {}'.format(i7_id, i7_seq))
                continue
            if suffix is None:
                # log.error('unknown fq suffix: {}'.format(fq))
                continue
            new_fq = os.path.join(i7_dir, name+suffix)
            if not file_exists(new_fq):
                symlink_file(fq, new_fq)
            # count reads
            fq_count = i7_fn_df.get(name, -1)
            if fq_count < 0:
                # fq_count = Fastx(new_fq).number_of_seq()
                if i7_seq not in self.i7_with_bc:
                    fq_count = Fastx(new_fq).number_of_seq()
                else:
                    fq_count = 1
                i7_fn_df.update({name:fq_count})
            if i7_seq not in self.i7_with_bc:
                if not suffix.startswith('_2'): # skip read2
                    log.info('check file: {} {}'.format(name, fq_count))
        # update read_count.json
        Config().dump(i7_fn_df, i7_fn)


    def demx_bc(self, x):
        """
        Parameters:
        ----------
        x : str
            barcode index table
        """
        args = self.__dict__.copy()
        i7_dir = os.path.join(self.out_dir, 'i7_index')
        i7_fq_list = list_fx(i7_dir, recursive=False)
        bc = file_prefix(x).replace('_table', '') # barcode seq
        bc_fq = [i for i in i7_fq_list if os.path.basename(i).startswith(bc)]
        args.update({
            'fq1': bc_fq[0] if len(bc_fq) > 0 else None,
            'fq2': bc_fq[1] if len(bc_fq) > 1 else None,
            'out_dir': os.path.join(self.out_dir, 'bc_index', bc),
            'index_type': 'barcode',
            'index_table': x,
        })
        DemxR1(**args).run() #


    def demx_i7_bc(self):
        args = self.__dict__.copy()
        tp, tb = self.split_i7_bc()
        if tp == 1:
            # 1 run i7
            self.rename_by_i7(tb)
        else:
            # 2. multi i7/bc mode
            tb_i7, tb_bc = tb #
            ## 2.1 run i7
            self.rename_by_i7(tb_i7)
            ## 2.2 run bc
            n_bc = 8 if len(tb_bc) > 8 else len(tb_bc)
            if n_bc > 0:
                with Pool(processes=n_bc) as pool:
                    pool.map(self.demx_bc, tb_bc)


    def wrap_i7_bc(self):
        tp, tb = self.split_i7_bc()
        # i7
        i7_df = Config().load(self.i7_fn_json) if file_exists(self.i7_fn_json) else None
        if i7_df is None:
            i7_df = {}
        i7_undemx = i7_df.get('undemx', 0)
        # barcode
        bc_df = {}
        i7_drop = []
        if tp > 1:
            for bc_f in tb[1]:
                bc = file_prefix(bc_f)
                bc = bc.replace('_table', '') # barcode seq
                bc_fc = os.path.join(self.out_dir, 'bc_index', bc, 'read_count.json')
                df = Config().load(bc_fc) #
                i7_undemx += df.get('undemx', 0)
                bc_df.update(df) # barcode files
                [i7_drop.append(k) for k,v in i7_df.items() if k.startswith(bc)]
        # merge
        i7_df.update(bc_df)
        [i7_df.pop(i, None) for i in i7_drop]
        i7_df.update({'undemx': i7_undemx}) # update undemx
        Config().dump(i7_df, self.fn_json)
        # i7+bc: rename files
        i7_dir = os.path.join(self.out_dir, 'i7_index')
        bc_dir = os.path.join(self.out_dir, 'bc_index')
        for i in list(self.idx.values()):
            # i7 index
            i7_list = list_fx2(i7_dir, i+'*', recursive=False)
            # bc index
            if file_exists(bc_dir):
                bc_list = list_fx2(bc_dir, i+'*', recursive=True)
            else:
                bc_list = []
            for a in i7_list + bc_list:
                a_new = os.path.join(self.out_dir, os.path.basename(a))
                if not file_exists(a_new):
                    symlink_file(a, a_new)


    def show_msg(self):
        msg = '\n'.join([
            '='*80,
            '{:>20} : {}'.format('Program', 'Demx_index'),
            '{:>20} : {}'.format('Date', get_date()),
            '{:>20} : {}'.format('fq_dir', self.data_dir),
            '{:>20} : {}'.format('index_table', self.index_table),
            '{:>20} : {}'.format('out_dir', self.out_dir),
            '{:>20} : {}'.format('report', self.report_txt),
            '{:>20} : {}'.format('mismatch', self.mismatch),
            '{:>20} : {}'.format('overwrite', 'yes' if self.overwrite else 'no'),
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
        # order
        i = 0
        f_stat = []
        # for k in list(self.idx.values()) + ['undemx']:
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
        self.show_msg()
        self.demx_i7_bc()
        self.wrap_i7_bc()
        self.report()


def get_args():
    example = '\n'.join([
        'De-multiplex single index/barcode',
        '1. only for i7 demx',
        '$ python demx2.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t i7',
        '2. only for barcode demx',
        '$ python demx2.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t barcode -x 2 -l 3 -r 2',
    ])
    parser = argparse.ArgumentParser(
        prog='demx2',
        description='De-multiplex single index',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--data_dir', dest='data_dir', required=True,
        help='Directory saving the fastq files')
    parser.add_argument('-s', '--index-table', dest='index_table',
        required=True,
        help='index table in csv format, [filename,barcode]')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    parser.add_argument('-x', '--barcode-in-read', dest='barcode_in_read',
        choices=[1, 2], default=2, type=int,
        help='barcode in read1/2, default: [2]')
    parser.add_argument('-l', '--barcode-n-left', type=int,
        dest='barcode_n_left', default=0,
        help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int,
        dest='barcode_n_right', default=0,
        help='bases locate on the right of barcode')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx2(**args).run()


if __name__ == '__main__':
    main()

#
