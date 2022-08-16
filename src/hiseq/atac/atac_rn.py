#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-3 (run _rn)
fastq files from: fq1/fq2
analysis-module:
"""


import os
# import sys
# import pathlib
# import shutil
from multiprocessing import Pool
from hiseq.utils.utils import log, update_obj, Config, init_cpu
from hiseq.utils.file import check_dir, symlink_file, file_abspath, fix_out_dir
from hiseq.utils.seq import check_fx_paired, fx_name
from hiseq.utils.hiseq_utils import list_hiseq_file
from hiseq.atac.atac_args import get_args_atac_rn
from hiseq.atac.atac_files import get_atac_dirs, get_atac_files
from hiseq.atac.atac_utils import (
    hiseq_call_peak, hiseq_bam2bw, hiseq_copy_r1, hiseq_merge_bam, hiseq_pcr_dup
)
from hiseq.atac.atac_qc import (
    qc_trim_summary, qc_align_summary, qc_lendist, qc_frip, 
    qc_tss_enrich, qc_genebody_enrich, qc_bam_cor, qc_peak_idr, 
    qc_peak_overlap, qc_bam_fingerprint
)
from hiseq.atac.atac_r1 import AtacR1
from hiseq.report.hiseq_report import HiSeqRpt


class AtacRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = AtacRnConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        # self.init_files()
        self.bam_list = list_hiseq_file(self.project_dir, 'bam_rmdup', 'r1')
        self.bw_list = list_hiseq_file(self.project_dir, 'bw', 'r1')
        self.peak_list = list_hiseq_file(self.project_dir, 'peak', 'r1')
        Config().dump(self.__dict__, self.config_yaml)


    def run_fx_r1(self, x):
        args = self.__dict__.copy()
        i = self.fq1.index(x) # index
        args.update({
            'fq1': x,
            'fq2': self.fq2[i],
            'smp_name': self.smp_name_list[i],
            'rep_list': None,
            'build_design': False,
            'design': None,
        })
        AtacR1(**args).run()


    def run_fx_rn(self):
        if self.parallel_jobs > 1 and len(self.fq1) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fx, self.fq1)
        else:
            [self.run_fx_r1(i) for i in self.fq1]


    def run(self):
        # 1. run r1
        self.run_fx_rn()
        # # 2. run rn
        # if len(self.rep_list) == 1:
        #     log.info('atac_rn() skipped, only 1 rep found')
        #     hiseq_copy_r1(self.project_dir)
        # else:   
        #     hiseq_merge_bam(self.project_dir, 'rn')
        #     hiseq_pcr_dup(self.project_dir, '_rn')
        #     hiseq_call_peak(self.project_dir, 'rn')
        #     hiseq_bam2bw(self.project_dir, 'rn')
        #     hiseq_call_peak(self.project_dir, 'rn')
        #     # qc_trim(self.project_dir, 'rn')
        #     # qc_align(self.project_dir, 'rn')
        #     qc_lendist(self.project_dir, 'rn')
        #     qc_frip(self.project_dir, 'rn')
        #     qc_tss_enrich(self.project_dir, 'rn')
        #     qc_genebody_enrich(self.project_dir, 'rn')
        #     # specific for rn
        #     qc_bam_cor(self.project_dir, 'rn')
        #     qc_peak_idr(self.project_dir, 'rn')
        #     qc_peak_overlap(self.project_dir, 'rn')
        #     qc_bam_fingerprint(self.project_dir, 'rn')
        # 4. report
        # HiSeqRpt(self.project_dir, overwrite=self.overwrite).run()


class AtacRnConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'fq1': None,
            'fq2': None,
            'out_dir': None,
            'smp_name': None,
            'smp_name_list': None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'atac_rn'
        self.out_dir = fix_out_dir(self.out_dir)
        self.init_fx()
        self.init_files()
        

    def init_fx(self):
        if not check_fx_paired(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not valid')
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        # fix smp_name (str, list)
        # if not isinstance(self.smp_name, str):
        self.smp_name = fx_name(self.fq1[0], fix_pe=True, fix_rep=True, fix_unmap=True)
        # for smp_name_list
        if isinstance(self.smp_name_list, list):
            stag = len(self.smp_name_list) == len(self.fq1)
        else:
            stag = False
        if not stag:
            self.smp_name_list = fx_name(self.fq1, fix_pe=True, fix_unmap=True)
        if self.smp_name == self.smp_name_list[0]:
            self.smp_name_list[0] += '_rep1'
        self.rep_list = [os.path.join(self.out_dir, i) for i in self.smp_name_list]


    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        atac_dirs = get_atac_dirs(self.out_dir, self.smp_name)
        atac_files = get_atac_files(self.out_dir, self.smp_name, self.fq1, self.fq2)
        self = update_obj(self, atac_dirs, force=True)
        self = update_obj(self, atac_files, force=True)
        # map(check_dir, atac_dirs.values())
        [check_dir(i) for i in atac_dirs.values()]


def get_args():
    return get_args_atac_rn()


def main():
    args = vars(get_args().parse_args())
    AtacRn(**args).run()


if __name__ == '__main__':
    main()

#