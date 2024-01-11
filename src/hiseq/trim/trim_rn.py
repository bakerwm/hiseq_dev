#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
required: cutadapt

Trim adapters, low quality bases, using cutadapt: multiple fx
1. Trim adapter,
2. remove dup
3. trim n-bases from read
3.1 TruSeq (NSR)
    - cut 9, -9; from both ends of read (subseq)
3.2 TruSeq (iCLIP)
    - cut 9; from read1
3.3 TruSeq (eCLIP)
    - cut 10; -7 from read1
    - cut 7, -10 from read2
optional
1. --rm-untrim, --save-too-short, ...
"""

import os

# import shutil
# import pathlib
import argparse
from multiprocessing import Pool
from hiseq.utils.seq import check_fx_args, check_fx_paired, fx_name
from hiseq.utils.file import check_dir, file_abspath
from hiseq.utils.utils import update_obj, Config
from hiseq.trim.cutadapt import get_args_trim
from hiseq.trim.trim_r1 import TrimR1, get_args_o


class TrimRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "fq1": None,
            "fq2": None,
            "smp_name": None,
            "parallel_jobs": 1,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "trim_rn"  #
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.init_fq()
        self.init_files()
        self.rep_list = [os.path.join(self.out_dir, i) for i in self.smp_name]
        Config().dump(self.__dict__.copy(), self.config_yaml)

    def init_fq(self):
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError("fq1, fq2 no valid")
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        if isinstance(self.smp_name, list):
            name_ok = len(self.fq1) == len(self.smp_name)
        else:
            name_ok = False
        if not name_ok:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)

    def init_files(self):
        self.project_dir = self.out_dir
        self.config_dir = os.path.join(self.project_dir, "config")
        self.config_yaml = os.path.join(self.config_dir, "config.yaml")
        check_dir(self.config_dir)

    def run_r1(self, x):
        args_local = self.__dict__.copy()
        i = self.fq1.index(x)  # index
        args_local.update(
            {
                "fq1": x,
                "fq2": self.fq2[i] if self.is_paired else None,
                # 'out_dir': self.out_dir,
                "smp_name": self.smp_name[i],
            }
        )
        TrimR1(**args_local).run()

    def run(self):
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_r1, self.fq1)
        else:
            [self.run_r1(fq) for fq in self.fq1]


def get_args_i():
    example = "\n".join(
        [
            "Examples:",
            "1. auto-detect adapters",
            "$ python trim_rn.py -1 fq1 -2 fq2 -o out_dir -m 20",
            "2. specific 3-apdapter",
            "$ python trim_rn.py -1 fq1 -2 fq2 -o out_dir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="trim",
        description="trim adapters",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        nargs="+",
        required=True,
        help="reads in FASTQ files, support (*.gz), 1 file.",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        nargs="+",
        default=None,
        help="The second file of pair-end reads",
    )
    return parser


def get_args():
    parser = get_args_o(get_args_i())  # I/O
    parser = get_args_trim(parser)
    parser.add_argument(
        "-j",
        "--parallel-jobs",
        type=int,
        dest="parallel_jobs",
        default=1,
        help="number of jobs run in parallel, default: [1]",
    )
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    # group.add_argument("-q", "--quiet", action="store_true")
    return parser


def main():
    args = vars(get_args().parse_args())
    TrimRn(**args).run()


if __name__ == "__main__":
    main()

#
