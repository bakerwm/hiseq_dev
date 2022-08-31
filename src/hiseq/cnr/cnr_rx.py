#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
CnR pipeline: level-2 (: rx; ip vs IgG)
loading fastq config from `design.yaml`, run pipeline, with specific parameters
analysis-module:
"""

import os
# import pathlib
# import argparse
from multiprocessing import Pool
from hiseq.utils.genome import Genome
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.utils.file import check_dir, file_abspath, file_exists, fix_out_dir, symlink_file
from hiseq.utils.seq import check_fx_args, fx_name
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu
from hiseq.utils.hiseq_utils import list_hiseq_file, read_hiseq
from hiseq.atac.atac_files import get_atac_dirs, get_atac_files
from hiseq.atac.atac_utils import (
    hiseq_merge_trim, hiseq_merge_bam, hiseq_copy_r1, hiseq_call_peak, 
    hiseq_bam2bw, hiseq_pcr_dup, hiseq_bw_compare
)
from hiseq.atac.atac_qc import (
    qc_trim_summary, qc_align_summary, qc_lendist, qc_frip, 
    qc_tss_enrich, qc_genebody_enrich, qc_bam_cor, qc_peak_idr, 
    qc_peak_overlap, qc_bam_fingerprint
)
from hiseq.cnr.cnr_args import get_args_cnr_rx
from hiseq.cnr.cnr_rn import CnrRn
from hiseq.report.hiseq_report import HiSeqRpt


class CnrRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = CnrRxConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)


    def run_smp_r1(self, is_ip=True):
        args = self.__dict__.copy()
        args.update({
            'build_design': False,
            'design': None,
            'fq_groups': None,
            'fq1': self.ip_fq1 if is_ip else self.input_fq1,
            'fq2': self.ip_fq2 if is_ip else self.input_fq2,
            'smp_name': self.ip_name if is_ip else self.input_name,
            'is_ip': is_ip,
            # 'parallel_jobs': 1,
        })
        CnrRn(**args).run()


    def run_smp_rn(self):
        args = self.__dict__.copy()
        self.run_smp_r1(is_ip=True)
        self.run_smp_r1(is_ip=False)


    def run_rx(self):
        """
        Run ip over input (IgG), for quality control
        """
        x = self.project_dir
        hiseq_call_peak(x, '_rx')
        # cnr_call_peak(x, 'rx') # call peak    
        if not self.fast_mode:
            hiseq_bw_compare(x, 'rx') # generate bw, ip.over.input    
            qc_tss_enrich(x, 'rx')
            qc_genebody_enrich(x, 'rx')
            qc_bam_cor(x, 'rx')
            qc_bam_fingerprint(x, hiseq_type='rx', bam_type='rn')


    def copy_ip(self):
        """
        copy ip files to rx directory: bam,bw,qc
        """
        # ip - bam,bw
        ip_dir = read_hiseq(self.ip_dir, 'rn')
        tx = [
            'bam', 'peak', 'peak_seacr', 'peak_seacr_top001', 'tss_enrich_png',
            'genebody_enrich_png', 'bam_fingerprint_png'
        ]
        for i in tx:
            symlink_file(list_hiseq_file(ip_dir, i), list_hiseq_file(self.project_dir, i))


    def prepare_files(self):
        # ip - bam,bw
        ra = read_hiseq(self.ip_dir, 'rn')
        symlink_file(ra.bam, self.ip_bam)
        symlink_file(ra.bw, self.ip_bw)
        # input - bam,bw
        rb = read_hiseq(self.input_dir, 'rn')
        if rb.is_hiseq:
            symlink_file(rb.bam, self.input_bam)
            symlink_file(rb.bw, self.input_bw)

            
    def run(self):
        # 1. run CnrRn
        self.run_smp_rn()
        # 2. run CnrRx
        self.prepare_files()
        if check_fx_args(self.input_fq1, self.input_fq2):
            self.run_rx()
        else:
            print('!ip - only, ...')
            self.copy_ip()
        # 3. generate report
        HiSeqRpt(self.project_dir, overwrite=self.overwrite).run()


class CnrRxConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
            'design': None,
            'build_design': False,
            'ip_fq1': None,
            'ip_fq2': None,
            'input_fq1': None,
            'input_fq2': None,
            'ip_name': None,
            'input_name': None,
            'smp_name': None,
            'out_dir': None,
            'parallel_jobs': 1,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_rx'
        self.out_dir = fix_out_dir(self.out_dir)
        self.init_index()
        self.init_name()
        self.init_files()
    
    
    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).fasize()
        else:
            raise ValueError('--genome or --extra-index; required')
        # update genome_size
        gs = 0
        with open(self.genome_size_file) as r:
            for line in r:
                gs += int(line.strip().split('\t')[1])
        self.genome_size = gs


    def init_name(self):
        """
        ip_name, input_name, smp_name
        ip_dir, input_dir
        smp_name: ip.vs.input
        """
        # ip_name: required
        if not isinstance(self.ip_name, str):
            self.ip_name = fx_name(self.ip_fq1[0], fix_pe=True, fix_rep=True)
        # input_name: optional, null
        if self.input_fq1 is None:
            self.input_name = 'null'
        else:
            if not isinstance(self.input_name, str):
                self.input_name = fx_name(self.input_fq1[0], 
                    fix_pe=True, fix_rep=True)
        # for ip.vs.input
        if not isinstance(self.smp_name, str):
            self.smp_name = '{}.vs.{}'.format(self.ip_name, self.input_name)


    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        atac_dirs = get_atac_dirs(self.out_dir, self.smp_name)
        atac_files = get_atac_files(self.out_dir, self.smp_name, self.fq1, self.fq2)
        self = update_obj(self, atac_dirs, force=True)
        self = update_obj(self, atac_files, force=True)
        cnr_files = {
            'ip_dir': os.path.join(self.out_dir, self.ip_name),
            'input_dir': os.path.join(self.out_dir, self.input_name),
            'ip_bam': os.path.join(self.bam_dir, self.ip_name+'.bam'),
            'input_bam': os.path.join(self.bam_dir, self.input_name+'.bam'),
            'ip_bw': os.path.join(self.bw_dir, self.ip_name+'.bigWig'),
            'input_bw': os.path.join(self.bw_dir, self.input_name+'.bigWig'),
        }
        self = update_obj(self, cnr_files, force=True)
        # map(check_dir, atac_dirs.values())
        [check_dir(i) for i in atac_dirs.values()]


def get_args():
    return get_args_cnr_rx()


def main():
    args = vars(get_args().parse_args())
    CnrRx(**args).run()


if __name__ == '__main__':
    main()

#