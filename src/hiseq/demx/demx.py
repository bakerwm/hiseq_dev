#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Demultiplex fastq by index and barcode; (double index)
# 1. split index table into subunits (i7:barcode)
# 2. if (i7) able to demx all: i7_only
# 3. if (i7) not able to demx; i7 + barcode
# position of i7/barcode:
i7: at the comment field of fastq record
barcode: at the beginning of read1/2
# alternative tools
1. defq, ultra fasta multi-threaded fq demultiplexing
https://github.com/OpenGene/defq
2. deML, maxlikelihood demultiplexing
https://github.com/grenaud/deML
# test
1: split by p7-index: index1_index2.r1.fq
2: split by inline barcode
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from collections import Counter
from contextlib import ExitStack
from hiseq.demx.demx_r1 import DemxR1
from hiseq.demx.sample_sheet import HiSeqIndex
from hiseq.utils.utils import log, update_obj, Config, get_date
from hiseq.utils.seq import list_fx, list_fx2, check_fx_args
from hiseq.utils.file import (
    check_dir, file_abspath, file_prefix, file_exists, symlink_file
)


class Demx(object):
    """
    Demultiplex fastq by i7+barcode
    Arguments
    ---------
    index_table : str
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
        'index_table': 'index.csv',
        'mismatch': 0,
        'demo': False,
        'gzipped': True,
    }
    >>> Demx2(**args).run()
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
            'fq1': None,
            'fq2': None,
            'index_type': 'i7', # i7, barcode
            'out_dir': None,
            'barcode_in_read': 2,
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
        self.hiseq_type = 'demx_r2'
        # fastq file
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq invalid: {}, {}'.format(self.fq1, self.fq2))
        # out_dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        if self.demo:
            self.out_dir = os.path.join(self.out_dir, 'demo')
        if not self.mismatch in range(4):
            raise ValueError('illegal mismatch: [{}], expect [0,1,2,3]'.format(
                self.mismatch))
        # index table
        self.init_files()
        self.idx = self.load_index(self.index_table, min_cols=4)
        # save config
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_files(self):
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.index_table = file_abspath(self.index_table)
        self.out_dir = file_abspath(self.out_dir)
        self.config_dir = os.path.join(self.out_dir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        self.fn_json = os.path.join(self.out_dir, 'read_count.json')
        self.report_txt = os.path.join(self.out_dir, 'report.txt')
        self.i7_table = os.path.join(self.out_dir, 'i7_table.csv')
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
                    px = [p1.match(i) is not None for i in idx_list]
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
                val = k1 if len(i7) == len(self.idx) else k2
                idx_list.append(','.join([v, val]))
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


    def demx_i7_bc(self):
        args = self.__dict__.copy()
        tp, tb = self.split_i7_bc()
        if tp == 1:
            # 1. single i7/bc mode
            args.update({
                'index_type': 'i7' if tb == self.i7_table else 'barcode',
                'index_table': tb,
            })
            DemxR1(**args).run() #
        else:
            # 2. multi i7/bc mode
            tb_i7, tb_bc = tb #
            ## 2.1 run i7
            i7_dir = os.path.join(self.out_dir, 'i7_index')
            args.update({
                'out_dir': i7_dir,
                'index_type': 'i7',
                'index_table': tb_i7,
            })
            DemxR1(**args).run() #
            ## 2.2 run bc
            i7_fq_list = list_fx(i7_dir, recursive=False)
            for bc_f in tb_bc:
                # print('!B-1', bc_f)
                bc = file_prefix(bc_f)
                bc = bc.replace('_table', '') # barcode seq
                bc_fq = [i for i in i7_fq_list if os.path.basename(i).startswith(bc)]
                # print('!A-9', bc_fq)
                args.update({
                    'fq1': bc_fq[0] if len(bc_fq) > 0 else None,
                    'fq2': bc_fq[1] if len(bc_fq) > 1 else None,
                    'out_dir': os.path.join(self.out_dir, 'bc_index', bc),
                    'index_type': 'barcode',
                    'index_table': bc_f,
                })
                DemxR1(**args).run() #


    def wrap_i7_bc(self):
        tp, tb = self.split_i7_bc()
        if tp == 1:
            return None # pass
        # i7+bc: read-count
        _, tb_bc = tb # list of barcode csv files
        # i7
        i7_fn = os.path.join(self.out_dir, 'i7_index', 'read_count.json')
        i7_df = Config().load(i7_fn)
        i7_undemx = i7_df.get('undemx', 0)
        # barcode
        bc_df = {}
        for bc_f in tb_bc:
            bc = file_prefix(bc_f)
            bc = bc.replace('_table', '') # barcode seq
            bc_fc = os.path.join(self.out_dir, 'bc_index', bc, 'read_count.json')
            df = Config().load(bc_fc) #
            i7_undemx += df.get('undemx', 0)
            bc_df.update(df) # barcode files
        # merge
        i7_df.update(bc_df)
        i7_df.update({'undemx': i7_undemx}) # update undemx
        Config().dump(i7_df, self.fn_json)
        # i7+bc: rename files
        i7_dir = os.path.join(self.out_dir, 'i7_index')
        bc_dir = os.path.join(self.out_dir, 'bc_index')
        for i in list(self.idx.values()):
            # i7 index
            i7_list = list_fx2(i7_dir, i+'*', recursive=False)
            [symlink_file(a, self.out_dir) for a in i7_list]
            # bc index
            bc_list = list_fx2(bc_dir, i+'*', recursive=True)
            [symlink_file(a, self.out_dir) for a in bc_list]


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
        # df = Config().load(self.fn_json) # index:count
        # total = sum(df.values())
        # if total < 1:
        #     total = 1000000 # default: 1M
        # # scale = 4e8/total # to_400M, deprecated
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
        self.demx_i7_bc()
        self.wrap_i7_bc()
        self.report()


def get_args():
    example = '\n'.join([
        'De-multiplex single index/barcode',
        '1. only for i7 demx',
        '$ python demx.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t i7',
        '2. only for barcode demx',
        '$ python demx.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t barcode -x 2 -l 3 -r 2',
    ])
    parser = argparse.ArgumentParser(
        prog='demx',
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
    parser.add_argument('-j', '--parallel-jobs', type=int, dest='parallel_jobs',
        default=1, help='number of jobs run in parallel, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demonstration, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx(**args).run()


if __name__ == '__main__':
    main()


# EOF
