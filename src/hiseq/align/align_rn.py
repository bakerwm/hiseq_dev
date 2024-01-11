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
from hiseq.utils.seq import Fastx, fx_name, check_fx_paired
from hiseq.utils.file import file_abspath, check_dir, list_dir, symlink_file
from hiseq.utils.hiseq_utils import list_hiseq_file, is_hiseq_dir
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.align.align_args import (
    get_args_io1,
    get_args_index2,
    get_args_align,
)
from hiseq.align.align_r1 import AlignR1
from hiseq.report.hiseq_report import HiSeqRpt


class AlignRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "aligner": "bowtie",
            "out_dir": None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "align_rn"
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.fx_format = Fastx(self.fq1).format  # fasta/q
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        if self.smp_name is None:
            self.smp_name = fx_name(
                self.fq1, fix_pe=self.is_paired, fix_unmap=True
            )
        self.init_index()
        self.init_files()
        Config().dump(self.__dict__.copy(), self.config_yaml)

    def init_index(self):
        self.index_list = check_index_args(**self.__dict__.copy())
        if len(self.index_list) == 0:
            raise ValueError("no index found")
        self.index_name = [AlignIndex(i).index_name() for i in self.index_list]
        # auto index_name, 01, 02, ...
        if len(self.index_name) > 1:
            self.index_name = [
                "{:02d}_{}".format(i, n)
                for i, n in zip(range(len(self.index_name)), self.index_name)
            ]

    def init_files(self):
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        self.config_dir = os.path.join(self.project_dir, "config")
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_args = {
            "config_yaml": os.path.join(self.config_dir, "config.yaml"),
            "bam": prefix + ".bam",
            "unmap": prefix + ".unmap." + self.fx_format,
            "unmap1": prefix + ".unmap.1." + self.fx_format,  #
            "unmap2": prefix + ".unmap.2." + self.fx_format,  #
            "align_stat": prefix + ".align.stat",
            "align_json": prefix + ".align.json",
            "align_flagstat": prefix + ".align.flagstat",
        }
        self = update_obj(self, default_args, force=True)
        check_dir(self.config_dir, create_dirs=True)

    def wrap_rn(self):
        """
        Save alignment to one file
        all files in project_dir
        """
        # save the last index as output
        last_dir = None
        for r1 in list_dir(self.project_dir, include_dirs=True):
            if is_hiseq_dir(r1, "r1"):
                last_dir = r1
        # save last index
        tx = [
            "bam",
            "align_stat",
            "align_json",
            "align_flagstat",
            "unmap",
            "unmap1",
            "unmap2",
        ]
        if is_hiseq_dir(last_dir, "r1"):
            for i in tx:
                src = list_hiseq_file(last_dir, i, "r1")
                dest = getattr(self, i)
                symlink_file(src, dest)

    def run(self):
        args = self.__dict__.copy()
        for index, index_name in zip(self.index_list, self.index_name):
            args.update(
                {
                    "index": index,
                    "index_name": index_name,
                    "index_list": None,
                    "rep_list": None,
                    "is_paired": None,
                }
            )
            bam, unmap1, unmap2 = AlignR1(**args).run()
            # update the unmap files for next round
            args.update(
                {
                    "fq1": unmap1,
                    "fq2": unmap2,
                }
            )
        # combine stat
        self.wrap_rn()
        # organize report
        HiSeqRpt(self.project_dir, overwrite=self.overwrite).run()


def get_args():
    parser = get_args_index2(get_args_io1())
    return get_args_align(parser)


def main():
    args = vars(get_args().parse_args())
    AlignRn(**args).run()


if __name__ == "__main__":
    main()
