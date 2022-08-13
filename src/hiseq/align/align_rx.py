#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
AlignR1() : fq=1, index=1, single fastq to single index
AlignRn() : fq=1, index=n, single fastq to multiple index
AlignRx() : fq=n, index=n, multiple fastq to multiple index
...
Working model:
Align {AlignConfig()}
  |- AlignRx()
  |- AlignRn()
  |- AlignR1()
AlignR1()
  |- {bowtie2, bowtie, STAR, bwa, ...}
"""

import os
import pathlib
import argparse
from hiseq.utils.utils import update_obj, Config
from hiseq.utils.seq import Fastx, fx_name, check_fx_paired, check_fx_args
from hiseq.utils.file import file_abspath, check_dir, list_dir, symlink_file
from hiseq.utils.hiseq_utils import list_hiseq_file, is_hiseq_dir
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.align.align_args import get_args_io2, get_args_index2, get_args_align
from hiseq.align.align_r1 import AlignR1
from hiseq.align.align_rn import AlignRn
from hiseq.report.hiseq_report import HiSeqRpt


class AlignRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie',
            'fq1': None,
            'fq2': None,
            'out_dir': None,
            'smp_name': None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'alignment_rx'
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.init_fq()
        self.init_files()  
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_fq(self):
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        if not check_fx_args(self.fq1, self.fq2, check_empty=True):
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # for names
        if isinstance(self.smp_name, list):
            stag = len(self.smp_name) == len(self.fq1)
        else:
            stag = False
        if not stag:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)
        self.rep_list = [os.path.join(self.out_dir, i) for i in self.smp_name]


    def init_files(self):
        self.project_dir = self.out_dir
        # self.project_dir = os.path.join(self.out_dir, self.smp_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        check_dir(self.config_dir)


    def run_fx_r1(self, x):
        args = self.__dict__.copy()
        i = self.fq1.index(x) # index
        args.update({
            'fq1': x,
            'fq2': self.fq2[i] if self.is_paired else None,
            'smp_name': self.smp_name[i],
            'rep_list': None,
            'is_paired': self.is_paired,
            })
        AlignRn(**args).run()


    def run(self):
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_r1, self.fq1)
        else:
            [self.run_r1(fq) for fq in self.fq1]


    def run(self):
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_fx_r1, self.fq1)
        else:
            [self.run_fx_r1(i) for i in self.fq1]
        # add project_dir
        self.project_dir = self.out_dir
        HiSeqRpt(self.project_dir, overwrite=self.overwrite).run()


def get_args():
    parser = get_args_index2(get_args_io2())
    return get_args_align(parser)


def main():
    args = vars(get_args().parse_args())
    AlignRx(**args).run()


if __name__ == '__main__':
    main()

#