#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-1 (main port)

mission-1: generate design.yaml

mission-2: run_pipe, parsing config from design.yaml
"""

import os
import pathlib
import argparse
from hiseq.atac.atac_rx import AtacRx
from hiseq.atac.atac_args import get_args_atac
from hiseq.utils.hiseq_utils import HiSeqDesignAtac
from hiseq.utils.utils import update_obj

# from hiseq.utils.file import file_abspath


class Atac(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)

    def run(self):
        if self.build_design:
            HiSeqDesignAtac(**self.__dict__).run()
        else:
            AtacRx(**self.__dict__).run()


def get_args():
    return get_args_atac()


def main():
    args = vars(get_args().parse_args())
    Atac(**args).run()


if __name__ == "__main__":
    main()

#
