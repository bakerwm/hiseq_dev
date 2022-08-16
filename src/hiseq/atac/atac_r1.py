#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
for single-pair fastq file, ATACseq analysis [group: 1, fx: 1]
analysis-module:
"""


import os
import sys
import pathlib
import shutil
from hiseq.utils.genome import Genome
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.utils.utils import log, update_obj, Config, init_cpu
from hiseq.utils.file import check_dir, symlink_file, file_abspath, fix_out_dir
from hiseq.utils.seq import check_fx_paired, fx_name
from hiseq.atac.atac_args import get_args_atac_r1 #, get_args_index1, get_args_atac1
from hiseq.atac.atac_files import get_atac_dirs, get_atac_files
from hiseq.atac.atac_utils import (
    hiseq_trim, hiseq_align_genome, hiseq_align_spikein, hiseq_call_peak, 
    hiseq_bam2bw, hiseq_pcr_dup
)
from hiseq.atac.atac_qc import (
    qc_trim_summary, qc_align_summary, qc_lendist, qc_frip,
    qc_tss_enrich, qc_genebody_enrich
)
from hiseq.report.hiseq_report import HiSeqRpt


class AtacR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args = AtacR1Config(**self.__dict__)
        self = update_obj(self, args.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)


    def run(self):
        x = self.project_dir
        # 1. symlink raw data
        raw_fq1, raw_fq2 = self.raw_fq_list
        symlink_file(self.fq1, raw_fq1, absolute_path=True)
        symlink_file(self.fq2, raw_fq2, absolute_path=True)
        # 2. run modules
        hiseq_trim(x, '_r1')
        hiseq_align_genome(x, '_r1')
        hiseq_pcr_dup(x, '_r1')
        # # sys.exit(1)
        hiseq_call_peak(x, '_r1')
        if isinstance(self.spikein_index, str):
            hiseq_align_spikein(x, '_r1')
        # 3. qc
        qc_trim_summary(x, '_r1')
        qc_align_summary(x, '_r1')
        qc_lendist(x, '_r1')
        qc_frip(x, '_r1')
        if not self.fast_mode:
            hiseq_bam2bw(x, '_r1') # optional 
            qc_tss_enrich(x, '_r1')
            qc_genebody_enrich(x, '_r1')
        # 4. report
        HiSeqRpt(x, overwrite=self.overwrite).run()


class AtacR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
            'fq1': None,
            'fq2': None,
            'out_dir': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'extra_index': None,
            'index_list': None,
            'smp_name': None,
            'rm_dup': True,
            'gene_bed': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'bin_size': 50,
            'normalizeUsing': 'RPGC',
            'genome_size': None,
            'genome_size_file': None,
            'genome_dir': None,
            'keep_tmp': None,
            'trimmed': False,
            'cut': False,
            'cut_to_length': 0,
            'recursive': False,
            'fast_mode': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'atac_r1'
        self.out_dir = fix_out_dir(self.out_dir)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        if self.gene_bed is None:
            self.gene_bed = Genome(self.genome).bed
        if self.cut:
            self.cut_to_length = 50
            self.recursive = True
        self.init_fq()
        self.init_index()
        self.init_files()
        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)


    def init_fq(self):
        if isinstance(self.fq1, str) and isinstance(self.fq2, str):
            if not check_fx_paired(self.fq1, self.fq2):
                raise ValueError('--fq1, --fq2 failed')
        else:
            raise ValueError('--fq1, --fq2, not str')
        # format
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = True # force PE reads
        # if not isinstance(self.smp_name, str):
        self.smp_name = fx_name(self.fq1, fix_pe=True, fix_unmap=True)
        self.rep_list = os.path.join(self.out_dir, self.smp_name)

  
    def init_index(self):
        index_list = check_index_args(**self.__dict__)
        if len(index_list) == 0:
            raise ValueError('no index found')
        # update: genome_size_file
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).fai
        else:
            raise ValueError('--genome or --extra-index; required')


    def init_files(self, create_dirs=True):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        atac_dirs = get_atac_dirs(self.out_dir, self.smp_name)
        atac_files = get_atac_files(self.out_dir, self.smp_name, self.fq1, self.fq2)
        self = update_obj(self, atac_dirs, force=True)
        self = update_obj(self, atac_files, force=True)
        # map(check_dir, atac_dirs.values())
        [check_dir(i) for i in atac_dirs.values()]


def get_args():
    return get_args_atac_r1()


def main():
    args = vars(get_args().parse_args())
    AtacR1(**args).run()


if __name__ == '__main__':
    main()

#