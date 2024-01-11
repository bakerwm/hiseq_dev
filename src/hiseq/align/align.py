#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
!!! switch from TOML to YAML, because of None type in python not saved in TOML format !!!
The structure of the classes:
Align()
AlignRx() : fq=n, index=n, multiple fastq to multiple index
AlignRn() : fq=1, index=n, single fastq to multiple index
AlignR1() : fq=1, index=1, single fastq to single index
AlignRp() : report()
AlignConfig()
AlignRxConfig()
AlignRnConfig()
AlignR1Config()
AlignRpConfig()
Working model:
Align {AlignConfig()}
  |- AlignRx()
  |- AlignRn()
  |- AlignR1()
AlignR1()
  |- {bowtie2, bowtie, STAR, bwa, ...}
"""

# import os
from hiseq.utils.utils import update_obj, get_date
from hiseq.align.align_rx import AlignRx, get_args


class Align(object):
    """
    The main port: fx: N; index: N
    support: multiple fx, multiple index
    force to list, even if `None`
    fq2: [None, None, ...]
    smp_name: [str, str, ...]
    index_name: [str, str, ...]
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "aligner": "bowtie2",
            "fq1": None,
            "fq2": None,
            "out_dir": None,
            "smp_name": None,
            "unique_only": False,
            "index_list": None,
            "extra_index": None,
            "genome": None,
            "genome_index": None,
            "spikein": None,
            "spikein_index": None,
            "to_rRNA": False,
            "rRNA_index": None,
            "to_chrM": False,
            "to_MT_trRNA": False,
            "threads": 1,
            "parallel_jobs": 1,
            "overwrite": False,
            "keep_tmp": False,
            "genome_size": 0,
            "genome_size_file": None,
            "extra_para": None,
            "verbose": False,
        }
        self = update_obj(self, args_init, force=False)

    def show_msg(self):
        msg = "\n".join(
            [
                "-" * 80,
                "{:>14s} : {}".format("Program", "Align"),
                "{:>14s} : {}".format("Date", get_date()),
                "{:>14s} : {}".format("aligner", self.aligner),
                "{:>14s} : {}".format("genome", self.genome),
                "{:>14s} : {}".format("index_list", self.index_list),
                "{:>14s} : {}".format("extra_index", self.extra_index),
                "{:>14s} : {}".format("fq1", self.fq1),
                "{:>14s} : {}".format("fq2", self.fq2),
                "{:>14s} : {}".format("out_dir", self.out_dir),
                "{:>14s} : {}".format("unique_only", self.unique_only),
                "{:>14s} : {}".format("extra_para", self.extra_para),
                "{:>14s} : {}".format("threads", self.threads),
                "{:>14s} : {}".format("parallel_jobs", self.parallel_jobs),
                "{:>14s} : {}".format("to_rRNA", self.to_rRNA),
                "{:>14s} : {}".format("to_chrM", self.to_chrM),
                "{:>14s} : {}".format("to_MT_trRNA", self.to_MT_trRNA),
                "-" * 80,
            ]
        )
        print(msg)

    def run(self):
        self.show_msg()
        AlignRx(**self.__dict__.copy()).run()


def main():
    args = vars(get_args().parse_args())
    Align(**args).run()


if __name__ == "__main__":
    main()


"""
The port for alignment
support: multi fx files, multi indexes

Align()
AlignRx() : fq=n, index=n, multiple fastq to multiple index
AlignRn() : fq=1, index=n, single fastq to multiple index
AlignR1() : fq=1, index=1, single fastq to single index
HiSeqRpt()

Working model:
Align {AlignConfig()}
  |- AlignRx()
  |- AlignRn()
  |- AlignR1()
  |- {bowtie2, bowtie, STAR, bwa, ...}

## requirements
- unique mapper
- multiple mapper
- config (yaml)
- saving unmapped ?!

## index
genome + {spikein, rRNA, MT}
extra_index {list}

top-level for alignment

input:
  - fq1: list (required)
  - fq2: list or None
  - genome: str or None
  - spikein: str or None
  - align-to-rRNA: bool
  - align-to-chrM: bool
  - smp_name: list or None (auto)
  - index_list: list or None
  - index_name: list or None (auto)
  - aligner: str (required)
  - index_list_equal: bool
  - n_map
  - unique_only
  - extra_para
  - parallel_jobs
  - threads

output:
  - bam: str
  - unmap1: str
  - unmap2: str or None 

priority:
  - 1. pickle
  - 2. genome, spikein, align-to-chrM, ...
  - x. extra_index

return bam, unmap1, unmap2

changelog
## update: 2021-03-08
1. run the codes independently
2. force 'bowtie2' to align rRNAs/chrM/...

## update: 2020-01-08
1. uniform code style: self ->  dict -> config.yaml

## update: 2020-04-24
1. rewrite the script, frame updated

to-do
1. force rRNA/chrM, unique + multiple
"""
