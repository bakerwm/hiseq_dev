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
from hiseq.utils.utils import update_obj
from hiseq.align.bowtie import Bowtie
from hiseq.align.bowtie2 import Bowtie2
from hiseq.align.hisat2 import Hisat2
from hiseq.align.star import Star
from hiseq.align.salmon import Salmon
from hiseq.utils.hiseq_utils import is_supported
from hiseq.align.align_index import AlignIndex, fetch_index
from hiseq.align.align_args import get_args_io1, get_args_index1, get_args_align


class AlignR1(object):
    """
    Alignment: fx: 1; index: 1
    return: bam, unmap-1, unmap-2
    save: log, stat, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_index()


    def init_index(self):
        args_init = {
            'aligner': 'bowtie',
            'index': None,
            'genome': None,
            'genome_path': None,
        }
        self = update_obj(self, args_init, force=False)
        # priority: index > genome
        if AlignIndex(self.index, self.aligner).is_valid():
            pass
        elif is_supported(self.genome, key='genome'):
            self.index = fetch_index(
                self.genome, aligner=self.aligner, genome_path=self.genome_path
            )
        else:
            raise ValueError('check -x, -g; no index found')


    def run(self):
        ad = {
            'bowtie': Bowtie,
            'bowtie2': Bowtie2,
            'star': Star,
            'salmon': Salmon,
            'hisat2': Hisat2,
#             'bwa': BWA,
#             'kallisto': Kallisto,
        }
        align = ad.get(self.aligner.lower(), None)
        if align is None:
            raise ValueError('unknown aligner: {}'.format(self.aligner))
        args = self.__dict__.copy()
        return align(**args).run() # bam, unmap-1, unmap-2
        # report ?


def get_args():
    parser = get_args_index1(get_args_io1())
    return get_args_align(parser)

def main():
    args = vars(get_args().parse_args())
    AlignR1(**args).run()


if __name__ == '__main__':
    main()
   