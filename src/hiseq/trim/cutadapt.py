#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run cutadapt v3.4 command
see https://cutadapt.readthedocs.io/ for documentation
Functions:
1. remove single adapter for SE reads (5', 3')
2. remove both adapter for PE reads (5', 3')
3. remove low quality bases on the end
4. remove extra bases
pass minus value to argument
"""


import os
import re
import shutil
import pathlib
import argparse
import pandas as pd
from hiseq.utils.seq import Fastx
from hiseq.utils.file import check_dir, file_exists, file_abspath
from hiseq.utils.seq import check_fx, check_fx_paired, fx_name
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd


class Cutadapt(object):
    """
    Wrap command for cutadapt
    Trimmer
    for cutadapt:
    1. adapters; quality; min_len
    1. trim adapter
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        self = update_obj(self, self.basic_args(), force=False)
        if not check_fx(self.fq1):
            raise ValueError('fq1 file not exists: {}'.format(self.fq1))
        self.hiseq_type = 'cutadapt_r1' #
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        self.init_fq()
        self.guess_adapter()
        self.prep_adapter()
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def basic_args(self):
        # cutadapt 3.4 with Python 3.8.6
        return {
            'library_type': None,
            'fq1': None,
            'fq2': None,
            'out_dir': None,
            'smp_name': None,
            # adapters
            'adapter3': None, # -a read1
            'adapter5': None, # -g read1
            'Adapter3': None, # -A read2
            'Adapter5': None, # -G read2
            'error_rate': 0.1, # -e
            'times': 1, # -n remove up to N adapters
            'overlap': 3, # -O overlap between read and adapter
            # modification
            'cut_before_trim': 0, # -u, --cut, Remove before trimming
            'qual_min': 20, # -q trim low-quality bases
            'cut_to_length': 0, # --length, shorten reads to length; after trimming
            'trim_n': False, # trim N's on the ends of reads
            # filtering
            'len_min': 15, # -m Discard reads shorter than LEN, 0
            'len_max': 0,  # -M Discard reads longer than Len, -1
            'max_n': 0.1, # --max-n, Discard reads with more than COUNT N bases
            'max_expected_errors': 0, # Discard reads exceed ERRORs
            'discard_trimmed': False, # --discard-trimmed
            'discard_untrimmed': False, # --discard-untrimmed
            # 'rm_untrim': False,
            # output
            # 'out_fq1': None, # -o first read in a pair
            # 'out_fq2 : None, # -p second read in a pair
            'save_too_short': False, # --too-short-output, --too-short-paired-output
            'save_too_long': False, # --too-long-output, --too-long-paired-output
            'save_untrimmed': False, # --untrimmed-output, --untrimmed-paired-output
            # pre-defined args
            'recursive': False, # trim adapters from begin of 3' adapter, end of 5' adapter
            'rm_polyN': False, # remove polyN, ['A', 'C', 'G', 'T', 'N', True]
            'threads': 4,
            'overwrite': False,
        }


    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        prefix = os.path.join(self.project_dir, self.project_name)
        # files
        default_files = {
            'config_yaml': os.path.join(self.config_dir,'config.yaml'),
            'cmd_sh': prefix + '.cmd.sh',
            'log': prefix + '.cutadapt.log',
            'trim_stat': prefix + '.cutadapt.trim.stat',
            'trim_yaml': prefix + '.cutadapt.trim.yaml',
            'trim_json': prefix + '.cutadapt.trim.json',
            'auto_adapter': prefix + '.auto_adapter.csv',
            'clean_fq': prefix + '.fq.gz', # update for SE
            'clean_fq1': prefix + '_1.fq.gz',
            'clean_fq2': prefix + '_2.fq.gz',
            'untrimmed_fq': prefix + '.untrimmed.fq.gz', # update for SE
            'untrimmed_fq1': prefix + '.untrimmed.1.fq.gz',
            'untrimmed_fq2': prefix + '.untrimmed.2.fq.gz',
            'too_short_fq': prefix + '.too_short.fq.gz', # update for SE
            'too_short_fq1': prefix + '.too_short.1.fq.gz',
            'too_short_fq2': prefix + '.too_short.2.fq.gz',
            'too_long_fq': prefix + '.too_long.fq.gz', # update for SE
            'too_long_fq1': prefix + '.too_long.1.fq.gz',
            'too_long_fq2': prefix + '.too_long.2.fq.gz',
        }
        self = update_obj(self, default_files, force=True) # key
        check_dir(self.config_dir)


    def init_fq(self):
        if not self.is_paired:
            # for SE, single-end
            # update filenames for se reads
            self.clean_fq1, self.untrimmed_fq1, self.too_short_fq1, self.too_long_fq1 = [
                self.clean_fq, self.untrimmed_fq, self.too_short_fq, self.too_long_fq
            ]
            self.clean_fq2, self.untrimmed_fq2, self.too_short_fq2, self.too_long_fq2 = [
                None, None, None, None
            ]


    def guess_adapter(self):
        """
        Guess the 3' adapters
        1. library_type
        2. from adapter3, Adapter3
        3. guess from the first 1M reads
        """
        lib = {
            'truseq': ['AGATCGGAAGAGC', 'AGATCGGAAGAGC'],
            'nextera': ['CTGTCTCTTATACA', 'CTGTCTCTTATACA'],
            'smallrna': ['TGGAATTCTCGGGTGCCAAGG', 'TGGAATTCTCGGGTGCCAAGG']
            }
        if self.library_type is None and self.adapter3 is None:
            log.info('Auto detect the adapters:')
            try:
                d = Fastx(self.fq1).detect_adapter()
                df = pd.DataFrame.from_dict(d, orient='index').sort_values(0,
                    ascending=False)
                self.library_type = df.index.to_list()[0] # first one
            except KeyError as e:
                log.error(e)
                self.library_type = 'truseq'
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
        elif isinstance(self.library_type, str):
            self.library_type = self.library_type.lower()
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
        else:
            pass


    def slice_str(self, x, step=1, min_size=10, direction=-1):
        """
        slice the string
        Parameters:
        -----------
        x : str
            the string
        step : int
            step size
        window : int
            width
        direction : int
            1 from left, -1 from right
        """
        if isinstance(x, str):
            x1 = x[::-1] if direction < 0 else x
            out = [x1[i:] for i in range(0, len(x1), step)]
            out = [i for i in out if len(i) >= min_size]
            return [i[::-1] for i in out] if direction < 0 else out


    def slice_adapter(self, x, recursive=False, **kwargs):
        step = kwargs.get('step', 1)
        min_size = kwargs.get('min_size', 8)
        direction = 1 # 1=from_left; -1=from_right
        if isinstance(x, str):
            if recursive:
                out = self.slice_str(x, step, min_size, direction)
            else:
                out = [x]
        elif isinstance(x, list):
            out = [j for i in x for j in self.slice_adapter(i, recursive, **kwargs)]
        else:
            out = []
        # filter ACTGN
        return [i for i in out if re.match('^[ACGTN]+$', i, re.IGNORECASE)] # list


    def prep_adapter(self):
        # to list; None => []
        # for 3' adapter
        self.adapter3 = self.slice_adapter(self.adapter3, self.recursive, direction=1)
        self.Adapter3 = self.slice_adapter(self.Adapter3, self.recursive, direction=1)
        # for 5' adapter
        self.adapter5 = self.slice_adapter(self.adapter5, self.recursive, direction=-1)
        self.Adapter5 = self.slice_adapter(self.Adapter5, self.recursive, direction=-1)


    def prep_cmd_adapter(self):
        cmd = []
        if len(self.adapter3) > 0:
            cmd += ['-a {}'.format(i) for i in self.adapter3]
        if len(self.Adapter3) > 0 and self.is_paired:
            cmd += ['-A {}'.format(i) for i in self.Adapter3]
        if len(self.adapter5) > 0:
            cmd += ['-g {}'.format(i) for i in self.adapter5]
        if len(self.Adapter5) > 0 and self.is_paired:
            cmd += ['-G {}'.format(i) for i in self.Adapter5]
        # remove polyN
        if isinstance(self.rm_polyN, str):
            cmd += ['-a {}{}'.format(self.rm_polyN, '{20}')]
        # matching
        cmd += [
            '-e {}'.format(self.error_rate),
            '-n {}'.format(self.times),
            '-O {}'.format(self.overlap),
        ]
        return ' '.join(cmd)


    def prep_cmd_additional(self):
        cmd = ['-q {} --trim-n'.format(self.qual_min)]
        if self.cut_before_trim:
            cmd += ['--cut {}'.format(self.cut_before_trim)]
        if self.cut_to_length:
            cmd += ['--length {}'.format(self.cut_to_length)]
        return ' '.join(cmd)


    def prep_cmd_filt(self):
        cmd = ['-m {}'.format(self.len_min)]
        if self.len_max > self.len_min:
            cmd += ['-M {}'.format(self.len_max)]
        # quality
        cmd += ['--max-n {}'.format(self.max_n)]
        if self.discard_untrimmed:
            cmd += ['--discard-untrimmed']
        elif self.discard_trimmed:
            cmd += ['--discard-trimmed']
        return ' '.join(cmd)


    def prep_cmd_output(self):
        cmd = []
        # short reads
        if self.is_paired:
            if self.save_too_short:
                cmd += [
                    '--too-short-output={}'.format(self.too_short_fq1),
                    '--too-short-paired-output={}'.format(self.too_short_fq2)
                ]
            if self.save_too_long:
                cmd += [
                    '--too-long-output={}'.format(self.too_long_fq1),
                    '--too-long-paired-output={}'.format(self.too_long_fq2)
                ]
            if self.save_untrimmed:
                cmd += [
                    '--untrimmed-output={}'.format(self.untrimmed_fq1),
                    '--untrimmed-paired-output={}'.format(self.untrimmed_fq2),
                ]
            cmd += ['-o {} -p {}'.format(self.clean_fq1, self.clean_fq2)]
        else:
            if self.save_too_short:
                cmd += ['--too-short-output={}'.format(self.too_short_fq)]
            if self.save_too_long:
                cmd += ['--too-long-output={}'.format(self.too_long_fq)]
            if self.save_untrimmed:
                cmd += ['--untrimmed-output={}'.format(self.untrimmed_fq)]
            cmd += ['-o {}'.format(self.clean_fq)]
        return ' '.join(cmd)


    def prep_cmd_io(self):
        cmd = ['-j {}'.format(self.threads)]
        if self.is_paired:
            cmd += ['{} {}'.format(self.fq1, self.fq2)]
        else:
            cmd += [self.fq1]
        cmd += ['1>{}'.format(self.log)]
        return ' '.join(cmd)


    def prep_cmd(self):
        cmd = ' '.join([
            'cutadapt',
            self.prep_cmd_adapter(),
            self.prep_cmd_additional(),
            self.prep_cmd_filt(),
            self.prep_cmd_output(),
            self.prep_cmd_io()
        ])
        return cmd


    def parse_log(self):
        """
        Parse the log file: .cutadapt.log
        SE:
        Total reads processed:               1,000,000
        Reads with adapters:                    78,522 (7.9%)
        Reads written (passing filters):       997,922 (99.8%)
        PE:
        Total read pairs processed:          1,000,000
          Read 1 with adapter:                  78,522 (7.9%)
          Read 2 with adapter:                  43,801 (4.4%)
        Pairs written (passing filters):       995,338 (99.5%)
        """
        p = re.compile('processed:\s+([\d,]+).*\n.*passing filters\):\s+([\d,]+)', re.DOTALL)
        n_total = 0
        n_clean = 0
        n_pct = 0
        if file_exists(self.log):
            with open(self.log, 'rt') as r:
                log_txt = r.read()
            m = p.search(log_txt) # match
            if m:
                n_total = int(m.group(1).replace(',', ''))
                n_clean = int(m.group(2).replace(',', ''))
                n_pct = '{:.2f}'.format(n_clean/n_total*100)
        header = ['#name', 'total', 'clean', 'percent']
        s = [self.smp_name, n_total, n_clean, n_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': n_total,
            'clean': n_clean,
            'percent': n_pct,
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # to global
        self.n_total = n_total
        self.n_clean = n_clean
        self.n_pct = n_pct
        return (n_total, n_clean, n_pct)


    def run(self):
        cmd = self.prep_cmd()
        with open(self.cmd_sh, 'wt') as w:
            w.write(cmd+'\n')
        # check output
        if self.is_paired:
            j = all(check_fx([self.clean_fq1, self.clean_fq2]))
        else:
            j = check_fx(self.clean_fq)
        if j and not self.overwrite:
            log.info('Cutadapt() skipped, file exists: {}'.format(self.out_dir))
        else:
            try:
                run_shell_cmd(cmd)
                self.parse_log()
            except:
                log.error('Cutadapt() failed, see: {}'.format(self.log))


def get_args_io():
    example = '\n'.join([
        'Examples:',
        '1. auto-detect adapters',
        '$ python cutadapt.py -1 fq1 -2 fq2 -o out_dir -m 20',
        '2. specific 3-apdapter',
        '$ python cutadapt.py -1 fq1 -2 fq2 -o out_dir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT',
    ])
    parser = argparse.ArgumentParser(
        prog='cutadapt',
        description='cutadapt: for single file',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
        help='FASTQ file, support (*.gz)')
    parser.add_argument('-2', '--fq2', default=None,
        help='The second of pair-end reads, or None')
    parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
        help='The directory to save results')
    parser.add_argument('-n', '--smp-name', dest='smp_name', default=None,
        help='The prefix of output files')
    return parser


def get_args_trim(parser):
    parser.add_argument('--library-type', dest='library_type', default=None,
        type=str, choices=['TruSeq', 'Nextera', 'smallRNA'],
        help='Type of the library structure, TruSeq, TruSeq standard library Nextera, Tn5 standard library, smallRNA, small RNA library')
    parser.add_argument('--rm-polyN', dest='rm_polyN',
        choices=['A', 'C', 'G', 'T'], default=None,
        help='remove polyN, default: [None]')
    parser.add_argument('-R', '--recursive', action='store_true',
        dest='recursive', help='trim adapter recursively')
    parser.add_argument('--times', type=int, default=1, dest='times',
        help='Remove up to COUNT adapters from read, default: [1]')
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', '--overlap', type=int, default=3,
        help='At least overlap between adapter and sequence, default: [3]')
    # add
    parser.add_argument('--cut-before-trim', dest='cut_before_trim', default='0',
        help='cut before trimming; example: 10: cut 10 from left; 0,-6; cut 6 from right; 9,-9, cut 9 from both ends, default [0]')
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')
    parser.add_argument('--cut-to-length', default=0, dest='cut_to_length',
        type=int,
        help='cut reads to from right, default: [0], full length')
    # filt
    parser.add_argument('-m', '--len-min', dest='len_min', default=15,
        type=int, help='Minimum length of reads after trimming, default [15]')
    parser.add_argument('-M', '--len-max', dest='len_max', default=0,
        type=int, help='Maximum length of reads after trimming, default [0], ignore')
    parser.add_argument('--rm-untrimmed', action='store_true',
        dest='discard_untrimmed',
        help='discard reads without adapter')
    parser.add_argument('--rm-trimmed', action='store_true',
        dest='discard_trimmed',
        help='discard reads with adapter')
    ## output
    parser.add_argument('--save-untrimmed', action='store_true',
        dest='save_untrimmed',
        help='Save untrimmed reads to file')
    parser.add_argument('--save-too-short', action='store_true',
        dest='save_too_short',
        help='Save too short reads to file')
    parser.add_argument('--save-too-long', action='store_true',
        dest='save_too_long',
        help='Save too short reads to file')
    # adapter
    parser.add_argument('-a', '--adapter3', default=None,
        help='3-Adapter sequence, default [None].')
    parser.add_argument('-g', '--adapter5', default=None,
        help='5-Adapter, default: None')
    parser.add_argument('-A', '--Adapter3', default=None,
        help='The 3 adapter of read2, default [None]')
    parser.add_argument('-G', '--Adapter5', default=None,
        help='The 5 adapter of read1, default: None')
    # additional
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='overwrite exists files')
    return parser


def main():
    args = vars(get_args_trim(get_args_io()))
    Cutadapt(**args).run()


if __name__ == '__main__':
    main()

#