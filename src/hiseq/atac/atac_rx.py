#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-2 (run pipe)

loading fastq config from `design.yaml`, run pipeline, with specific parameters

analysis-module:
"""

import os
# import pathlib
# import argparse
from multiprocessing import Pool
from hiseq.utils.genome import Genome
from hiseq.utils.file import check_dir, file_abspath, file_exists, fix_out_dir
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu
from hiseq.atac.atac_rn import AtacRn
from hiseq.report.hiseq_report import HiSeqRpt


class AtacRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args = AtacRxConfig(**self.__dict__)
        self = update_obj(self, args.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)


    def run_group_r1(self, x):
        dx = self.fq_groups.get(x, {})
        args = self.__dict__.copy()
        args.update({
            'build_design': False,
            'design': None,
            'fq_groups': None,
            'fq1': [v[0] for k,v in dx.items()],
            'fq2': [v[1] for k,v in dx.items()],
        })
        if len(self.fq_groups) > 1:
            args['parallel_jobs'] = 1 # force
        AtacRn(**args).run()


    def run_group_rx(self):
        if self.parallel_jobs > 1 and len(self.fq_groups) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_group_r1, self.fq_groups)
        else:
            [self.run_group_r1(i) for i in self.fq_groups]


    def run(self):
        # 1. run AtacRn->AtacR1
        self.run_group_rx()
        # 2. report
        HiSeqRpt(self.project_dir, overwrite=self.overwrite).run()


class AtacRxConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'design': None, # required
            'build_design': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'atac_rx'
        self.out_dir = fix_out_dir(self.out_dir)
        if not file_exists(self.design):
            raise ValueError('file not exists, --design {}'.format(self.design))
        self.init_files()
        self.init_fx()


    def init_files(self):
        # dirs
        self.project_dir = os.path.join(self.out_dir)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.report_dir = os.path.join(self.project_dir, 'report')
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
        }
        self = update_obj(self, default_files, force=True) # key
        check_dir([self.config_dir, self.report_dir])


    def init_fx(self):
        """
        Loading fastq config from design
        """
        self.fq_groups = Config().load(self.design)
        if len(self.fq_groups) == 0:
            raise ValueError('no fq found: {}'.format(self.design))


def get_args():
    return get_args_atac_rx()


def main():
    args = vars(get_args().parse_args())
    AtacRx(**args).run()


if __name__ == '__main__':
    main()


#