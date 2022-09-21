#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
for single-pair fastq file, CUT&Tag, CUT&Run analysis [group: 1, fx: 1]
fq: PE or SE
    - ip (required)
    - input/IgG (optional)
1. trimming
    double-trimming, recursive
    a. full length adapter 
    b. <6bp adapters
2. Alignment
    dovetail alignment using bowtie2
    --dovetail
    --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700
3. pcr-dup: 
    - remove 
    - keep
4. peak calling
    - SEACR
    - Macs2
5. cut metrix
6. enrichment of motif 
"""


import os
import sys
import pathlib
import shutil
from hiseq.utils.genome import Genome
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.utils.utils import log, update_obj, Config, init_cpu
from hiseq.utils.file import check_dir, symlink_file, file_abspath, fix_out_dir
from hiseq.utils.seq import check_fx_paired, fx_name, check_fx_args
from hiseq.cnr.cnr_args import get_args_fast
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


class CnrR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args = CnrR1Config(**self.__dict__)
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


class CnrR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
            'is_ip': True,
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
        self.hiseq_type = 'cnr_r1'
        self.out_dir = fix_out_dir(self.out_dir)
        if self.gene_bed is None:
            self.gene_bed = Genome(self.genome).bed
        if self.cut:
            self.cut_to_length = 50
            self.recursive = True
        self.init_fx()
        self.init_index()
        self.init_files()
        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)


    def init_fx(self):
        """
        1. PE or SE
        """
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not valid, check above message')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)


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
        # update genome_size
        gs = 0
        with open(self.genome_size_file) as r:
            for line in r:
                gs += int(line.strip().split('\t')[1])
        self.genome_size = gs


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
    parser = get_args_atac_r1()
    return get_args_fast(parser)


def main():
    args = vars(get_args().parse_args())
    CnrR1(**args).run()


if __name__ == '__main__':
    main()

#