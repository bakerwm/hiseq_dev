#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sample sequences by number; -n
Generate sample file from fastx file
randome: seqkit sample
head:
tail:
"""


import os
import re
import glob
import argparse
from multiprocessing import Pool
from hiseq.utils.seq import Fastx
from hiseq.utils.utils import update_obj, log
from hiseq.utils.file import file_exists, file_abspath, check_dir


class Sample(object):
    """
    Get subset of fastx file
    random: true|false
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'fx': None,
            'number': None,
            'out_dir': None,
            'random': False,
            'overwrite': False,
            'parallel_jobs': 1 
        }
        self = update_obj(self, args_init, force=False)
        # file exists
        if isinstance(self.fx, str):
            self.fx = [self.fx]
        elif isinstance(self.fx, list):
            pass
        else:
            raise ValueError('-i expect str or list, got {}'.format(
                type(self.fx).__name__))
        # number
        if isinstance(self.number, int):
            self.number = abs(self.number)
        else:
            raise ValueError('-n expect int, got {}'.format(
                type(self.number).__name__))
        # output
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        check_dir(self.out_dir, create_dirs=True)


    def sampleR1(self, fx):
        fx_out = os.path.join(self.out_dir, os.path.basename(fx))
        if not fx_out.endswith('.gz'):
            fx_out += '.gz' # gzip output
        if file_exists(fx_out) and not self.overwrite:
            log.info('sample() skipped, file exists: {}'.format(fx_out))
        else:
            if self.random:
                log.info('extract {} random records from: {}'.format(self.number, fx))
                Fastx(fx).sample_random(fx_out, n=self.number)
            else:
                log.info('extract {} records from: {}'.format(self.number, fx))
                Fastx(fx).sample(fx_out, n=self.number)


    def sampleRn(self):
        # run each fq
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.sampleR1, self.fx)


    def run(self):
        self.sampleRn()


def get_args():
    parser = argparse.ArgumentParser(description='hiseq sample')
    parser.add_argument('-i', '--fx', nargs='+', required=True,
        help='fastx files')
    parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
        help='output directory to save results')
    parser.add_argument('-n', '--number', type=int, default=1000,
        help='Number of records, default: 1000')
    parser.add_argument('-r', '--random', action='store_true',
        help='Get random subset records'),
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite the exists files')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1, type=int,
        help='Number of threads run in parallel, default [1]')
    return parser


def main():
    args = vars(get_args().parse_args())
    Sample(**args).run()


if __name__ == '__main__':
    main()

#