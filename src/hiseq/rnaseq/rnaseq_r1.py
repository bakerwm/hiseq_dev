#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
RNAseq for single fq, single index: r1
"""

import os
import sys
import pathlib
import argparse
# from hiseq.rnaseq.rnaseq_rp import RnaseqRp
# from hiseq.rnaseq.utils import (
#     hiseq_trim, hiseq_align_spikein, rnaseq_align_rRNA, hiseq_align_genome,
#     rnaseq_quant, hiseq_bam2bw, qc_trim_summary, qc_align_summary, 
#     qc_genebody_enrich
# )
# from hiseq.utils.file import (
#     check_dir, check_fx_paired, symlink_file, file_exists, file_abspath, 
#     file_prefix, fx_name, check_fx_args
# )
# from hiseq.utils.utils import (
#     log, update_obj, Config, get_date, init_cpu, read_hiseq, is_supported, 
#     print_dict
# )

####
from hiseq.utils.utils import log, update_obj, Config, init_cpu
from hiseq.utils.file import check_dir, symlink_file, file_abspath, fix_out_dir
from hiseq.utils.seq import check_fx_paired, fx_name, check_fx_args
from hiseq.align.align_index import AlignIndex, check_index_args, fetch_index
from hiseq.utils.genome import Genome
from hiseq.report.hiseq_report import HiSeqRpt
from hiseq.atac.atac_utils import (
    hiseq_trim, hiseq_align_genome, hiseq_align_spikein, hiseq_bam2bw
)
from hiseq.atac.atac_qc import (
    qc_trim_summary, qc_align_summary, qc_tss_enrich, qc_genebody_enrich
)
from hiseq.rnaseq.utils import rnaseq_quant, rnaseq_align_rRNA
from hiseq.rnaseq.rnaseq_args import get_args_rnaseq_r1, get_args_fast
from hiseq.rnaseq.rnaseq_files import get_rnaseq_dirs, get_rnaseq_files


class RnaseqR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = RnaseqR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)

        
    def run(self):
        x = self.project_dir
        # 1. symlink raw data
        raw_fq1, raw_fq2 = self.raw_fq_list
        symlink_file(self.fq1, raw_fq1, absolute_path=True)
        symlink_file(self.fq2, raw_fq2, absolute_path=True)
        # 2. run modules
        hiseq_trim(x, '_r1')
        hiseq_align_spikein(x, '_r1')
        rnaseq_align_rRNA(x, 'r1')
        hiseq_align_genome(x, 'r1')
        rnaseq_quant(x, 'r1')
        # 3. qc
        qc_trim_summary(self.project_dir, '_r1')
        qc_align_summary(self.project_dir, '_r1')
        if not self.fast_mode:
            hiseq_bam2bw(x, '_r1') # optional 
            qc_tss_enrich(x, '_r1')
            qc_genebody_enrich(x, '_r1')
        # 4. report
        HiSeqRpt(x, overwrite=self.overwrite).run()


class RnaseqR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'is_mut': False,
            'aligner': 'STAR',            
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
            'to_rRNA': False,
            'rRNA_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'verbose': False,
            'bin_size': 50,
            'genome_size': 0,
            'genome_size_file': None,
            'genome_dir': None,
            'gene_bed': None,
            'keep_tmp': False,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False,
            'extra_para': None,
            'norm_project': None,
            'fast_mode': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_r1'
        self.out_dir = fix_out_dir(self.out_dir)
        if self.gene_bed is None:
            self.gene_bed = Genome(self.genome).bed
        if self.cut:
            self.cut_to_length = 50
            self.recursive = True
        # # output
        # if self.out_dir is None:
        #     self.out_dir = str(pathlib.Path.cwd())
        # self.out_dir = file_abspath(self.out_dir)
        # if self.gene_bed is None:
        #     self.gene_bed = Genome(self.genome).gene_bed('ensembl')
        # if self.gene_gtf is None:
        #     self.gene_gtf = Genome(self.genome).gene_gtf('ensembl')
        self.init_fx()
        self.init_index()
        self.init_files()
        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)        
        # if not isinstance(self.smp_name, str):
        #     self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)
        # self.init_files()
        # self.init_index()
        # self.threads, _ = init_cpu(self.threads, 1)
        # Config().dump(self.__dict__, self.config_yaml)


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


    # update: genome_size_file
    def init_index(self):
        index_list = check_index_args(**self.__dict__)
        if len(index_list) == 0:
            raise ValueError('no index found')
        # update: genome_size_file          
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


    def init_files(self, create_dirs=True):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        atac_dirs = get_rnaseq_dirs(self.out_dir, self.smp_name)
        atac_files = get_rnaseq_files(self.out_dir, self.smp_name, self.fq1, self.fq2)
        self = update_obj(self, atac_dirs, force=True)
        self = update_obj(self, atac_files, force=True)
        # map(check_dir, atac_dirs.values())
        [check_dir(i) for i in atac_dirs.values()]


def get_args():
    parser = get_args_rnaseq_r1()
    return get_args_fast(parser)

    # def init_files(self):
    #     # dirs
    #     self.project_name = self.smp_name
    #     self.project_dir = os.path.join(self.out_dir, self.project_name)
    #     self.config_dir = os.path.join(self.project_dir, 'config')
    #     default_dirs = {
    #         'raw_dir': 'raw_data',
    #         'clean_dir': 'clean_data',
    #         'align_dir': 'align',
    #         'spikein_dir': 'spikein',
    #         'rRNA_dir': 'rRNA',
    #         'bam_dir': 'bam_files',
    #         'bg_dir': 'bg_files',
    #         'bw_dir': 'bw_files',
    #         'count_dir': 'count',
    #         'qc_dir': 'qc',
    #         'report_dir': 'report'
    #     }
    #     for k, v in default_dirs.items():
    #         default_dirs[k] = os.path.join(self.project_dir, v)
    #     self = update_obj(self, default_dirs, force=True) # key
    #     # files        
    #     trim_prefix = os.path.join(self.clean_dir, self.smp_name)
    #     spikein_prefix = os.path.join(self.spikein_dir, self.smp_name)
    #     rRNA_prefix = os.path.join(self.rRNA_dir, self.smp_name)
    #     align_prefix = os.path.join(self.align_dir, self.smp_name)
    #     count_prefix = os.path.join(self.count_dir, self.smp_name)
    #     default_files = {
    #         # basic files
    #         'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
    #         'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
    #         'bam': os.path.join(self.bam_dir, self.smp_name+'.bam'),
    #         'bw': os.path.join(self.bw_dir, self.smp_name+'.bigWig'),
    #         'bw_fwd': os.path.join(self.bw_dir, self.smp_name+'.fwd.bigWig'),
    #         'bw_rev': os.path.join(self.bw_dir, self.smp_name+'.rev.bigWig'),
            
    #         # trimming
    #         'trim_stat': trim_prefix+'.trim.stat',
    #         'trim_json': trim_prefix+'.trim.json',
            
    #         # quantification files
    #         'strandness_json': count_prefix+'.strandness.json',
    #         'count_sens': count_prefix+'.sens.txt',
    #         'count_anti': count_prefix+'.anti.txt',           
            
    #         # align files (genome)
    #         'align_scale_json': self.bam_dir + '/' + 'scale.json',
    #         # 'align_scale_json': align_prefix+'.scale.json',
    #         'align_stat': align_prefix+'.align.stat',
    #         'align_json': align_prefix+'.align.json',
    #         'align_flagstat': align_prefix+'.align.flagstat',
    #         'unmap': align_prefix+'.unmap.fastq',
    #         'unmap1': align_prefix+'.unmap.1.fastq',
    #         'unmap2': align_prefix+'.unmap.2.fastq',
            
    #         # spikein files
    #         'spikein_bam': spikein_prefix+'.bam',
    #         'spikein_scale_json': spikein_prefix+'.scale.json',
    #         'spikein_stat': spikein_prefix+'.align.stat',
    #         'spikein_json': spikein_prefix+'.align.json',
    #         'spikein_flagstat': spikein_prefix+'.align.flagstat',
    #         'spikein_unmap': spikein_prefix+'.unmap.fastq',
    #         'spikein_unmap1': spikein_prefix+'.unmap.1.fastq',
    #         'spikein_unmap2': spikein_prefix+'.unmap.2.fastq',
            
    #         # rRNA files
    #         'rRNA_bam': rRNA_prefix+'.bam',
    #         'rRNA_stat': rRNA_prefix+'.align.stat',
    #         'rRNA_json': rRNA_prefix+'.align.json',
    #         'rRNA_flagstat': rRNA_prefix+'.align.flagstat',
    #         'rRNA_unmap': rRNA_prefix+'.unmap.fastq',
    #         'rRNA_unmap1': rRNA_prefix+'.unmap.1.fastq',
    #         'rRNA_unmap2': rRNA_prefix+'.unmap.2.fastq',
            
    #         # qc
    #         'trim_summary_json': os.path.join(self.qc_dir, '00.trim_summary.json'),
    #         'align_summary_json': os.path.join(self.qc_dir, '01.alignment_summary.json'),
    #         'genebody_enrich_matrix': os.path.join(self.qc_dir, '05.genebody_enrich.mat.gz'),
    #         'genebody_enrich_matrix_log': os.path.join(self.qc_dir, '05.genebody_enrich.log'),
    #         'genebody_enrich_png': os.path.join(self.qc_dir, '05.genebody_enrich.png'),
    #         'genebody_enrich_cmd': os.path.join(self.qc_dir, '05.genebody_enrich.cmd.sh'),
    #         'bam_cor_npz': os.path.join(self.qc_dir, '06.bam_cor.npz'),
    #         'bam_cor_counts': os.path.join(self.qc_dir, '06.bam_cor.counts.tab'),
    #         'bam_cor_heatmap_png': os.path.join(self.qc_dir, '06.bam_cor.cor_heatmap.png'),
    #         'bam_cor_pca_png': os.path.join(self.qc_dir, '06.bam_cor.cor_PCA.png')
    #     }
    #     self = update_obj(self, default_files, force=True) # key
    #     # update fq1,fq2 files: support for SE
    #     if self.is_paired:
    #         self.raw_fq_list = [
    #             self.raw_dir + '/' + os.path.basename(self.fq1),
    #             self.raw_dir + '/' + os.path.basename(self.fq2),
    #         ]
    #         self.clean_fq_list = [
    #             self.clean_dir + '/' + os.path.basename(self.fq1),
    #             self.clean_dir + '/' + os.path.basename(self.fq2),
    #         ]
    #     else:
    #         self.raw_fq_list = [
    #             self.raw_dir + '/' + os.path.basename(self.fq1),
    #             None,
    #         ]
    #         self.clean_fq_list = [
    #             self.clean_dir + '/' + os.path.basename(self.fq1),
    #             None,
    #         ]
    #         self.spikein_unmap1, self.spikein_unmap2 = [self.spikein_unmap, None]
    #         self.rRNA_unmap1, self.rRNA_unmap2 = [self.rRNA_unmap, None]
    #         self.unmap1, self.unmap2 = [self.unmap, None]
    #     # create dirs
    #     dir_list = [
    #         self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
    #         self.spikein_dir, self.rRNA_dir, self.bam_dir, self.bg_dir,
    #         self.bw_dir, self.count_dir, self.qc_dir, self.report_dir
    #     ]
    #     check_dir(dir_list)



        
def get_args():
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python rnaseq_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6',
        '2. for specific index',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te',
    ])    
    parser = argparse.ArgumentParser(
        prog='rnaseq_r1',
        description='run rnaseq_r1',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = argparse.ArgumentParser(
        description='RNA-seq pipeline')
    parser.add_argument('-1', '--fq1', default=None, required=True,
        help='read1 files, (or read1 of PE reads)')
    parser.add_argument('-2', '--fq2', default=None,
        help='read2 of PE reads')
    parser.add_argument('-o', '--out_dir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: hg38')
    parser.add_argument('--gtf', dest='gene_gtf', default=None,
        help='The gtf file for quantification, defaut: genome.gtf (None)')
    parser.add_argument('--bed', dest='gene_bed', default=None,
        help='The BED of genes')
    # optional arguments - 0
    parser.add_argument('-n', '--smp-name', dest='smp_name', default=None,
        help='Name of the samples, default: from fq1 prefix')
    # optional arguments - 0
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    # optional arguments - 1
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    ## extra: index
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--genome-index', dest="genome_index", default=None,
        help='align index of genome')
    parser.add_argument('-k', '--spikein', default=None,
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--spikein-index', dest="spikein_index", default=None,
        help='align index of spikein')    
    parser.add_argument('--to-rRNA', dest='to_rRNA', action='store_true',
        help='Align to rRNA')
    parser.add_argument('--rRNA-index', dest="rRNA_index", default=None,
        help='align index of rRNA')   
    parser.add_argument('--aligner', default='STAR',
        choices=['STAR', 'bowtie', 'bowtie2', 'bwa', 'hisat2', 'kallisto',
            'salmon'],
        help='Aligner option: [STAR, bowtie, bowtie2, bwa], default: [STAR]')
    ## extra:
    parser.add_argument('-bs', '--bin-size', default=10, type=int,
        help='bin size of the bigWig file, default [10]')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. \
        default: [None]')
    parser.add_argument('--norm-project', dest='norm_project', default=None,
        help='The RNAseq_Rx project, for parseing norm scale. eg: \
        RNAseq_gene/wt.vs.mut for RNAseq_te, default: [None]')
    return parser


def main():
    args = vars(get_args().parse_args())
    RnaseqR1(**args).run()


if __name__ == '__main__':
    main()

#
