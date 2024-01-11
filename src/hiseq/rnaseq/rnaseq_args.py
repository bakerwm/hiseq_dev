#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import re
import argparse
from hiseq.rnaseq.rnaseq_args import get_args_index1, get_args_opt


def get_args_cnr():
    parser = get_args_io4()
    # parser = get_args_fast(parser)
    parser = get_args_index1(parser)
    return get_args_opt(parser)


def get_args_cnr_rx():
    parser = get_args_io3()
    # parser = get_args_fast(parser)
    parser = get_args_index1(parser)
    return get_args_opt(parser)


def get_args_fast(parser):
    parser.add_argument(
        "--fast",
        dest="fast_mode",
        action="store_true",
        help="run in fast mode, for a quick look of the data",
    )
    return parser


def get_args_atac_r1():
    parser = get_args_io1()
    parser = get_args_index1(parser)
    return get_args_opt(parser)


def get_args_io1():
    example = "\n".join(
        [
            "Examples:",
            "1. support fastq input",
            "$ python rnaseq_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6",
            "2. for specific index",
            "$ python rnaseq_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="rnaseq_r1",
        description="rnaseq_r1: for single PE reads",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1", "--fq1", required=True, help="Fastq file, read1 of PE"
    )
    parser.add_argument(
        "-2", "--fq2", required=True, help="Fastq file, read2 of PE"
    )
    # parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
    #     help='Directory saving results, default: [pwd]')
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="The name of the sample",
    )
    return parser


def get_args_io2():
    example = "\n".join(
        [
            "Examples:",
            "1. support fastq input",
            "$ python rnaseq_rn.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6",
            "2. for specific index",
            "$ python rnaseq_rn.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="rnaseq_rn",
        description="rnaseq_r1: for single PE reads",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1", "--fq1", nargs="+", required=True, help="Fastq file, read1 of PE"
    )
    parser.add_argument(
        "-2", "--fq2", nargs="+", required=True, help="Fastq file, read2 of PE"
    )
    # parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
    #     help='Directory saving results, default: [pwd]')
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="The name of the sample",
    )
    parser.add_argument(
        "-N",
        "--smp-name-list",
        dest="smp_name_list",
        default=None,
        help="The list of names for fastq files",
    )
    return parser


def get_args_io3():
    example = "\n".join(
        [
            "Examples:",
            "1. support fastq input",
            "$ python rnaseq_rx.py -o results -g dm6 --wt-fq1 --wt-fq2 --mut-fq1 --mut-fq2",
            "2. for specific index",
            "$ python rnaseq_rx.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="rnaseq_r1",
        description="run rnaseq",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser = argparse.ArgumentParser(description="RNA-seq pipeline")
    parser.add_argument(
        "-m1",
        "--mut-fq1",
        nargs="+",
        required=True,
        dest="mut_fq1",
        help="read1 files, (or read1 of PE reads) of treatment/mutant samples",
    )
    parser.add_argument(
        "-m2",
        "--mut-fq2",
        nargs="+",
        required=False,
        dest="mut_fq2",
        help="read2 files, (or read2 of PE reads) of treatment/mutant samples",
    )
    parser.add_argument(
        "-w1",
        "--wt-fq1",
        nargs="+",
        required=True,
        dest="wt_fq1",
        help="read1 files, (or read1 of PE reads) of control/wt samples",
    )
    parser.add_argument(
        "-w2",
        "--wt-fq2",
        nargs="+",
        required=False,
        dest="wt_fq2",
        help="read2 files, (or read2 of PE reads) of control/wt samples",
    )
    # optional arguments - 0
    parser.add_argument(
        "-mn",
        "--mut-name",
        dest="mut_name",
        default=None,
        help="Name of mutant samples",
    )
    parser.add_argument(
        "-wn",
        "--wt-name",
        dest="wt_name",
        default=None,
        help="Name of wildtype/control samples",
    )
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="Name of the samples, default: from fq1 prefix",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default=None,
        dest="out_dir",
        help="The directory to save results, default, \
        current working directory.",
    )
    return parser


def get_args_io4():
    example = "\n".join(
        [
            "Examples:",
            "1. support fastq input",
            "$ python rnase.py -b -d rnaseq.json -r fq_dir --mut dTAG --wt DMSO",
            "2. for specific index",
            "$ python rnaseq.py -d rnaseq.json -o results -g dm6",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="cnr",
        description="cnr: for multiple groups (ip vs input)",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-b",
        "--build-design",
        dest="build_design",
        action="store_true",
        help="generate design.yaml, with --fq-dir, or fq1/fq2",
    )
    parser.add_argument(
        "-d",
        "--design",
        required=True,
        help="The file saving fastq files config; generated by cnr_rd.py",
    )
    parser.add_argument(
        "-r",
        "--fq-dir",
        dest="fq_dir",
        help="Path to the directory, contains fastq files, eg: _rep1_1.fq.gz",
    )
    parser.add_argument(
        "-1",
        "--fq1",
        nargs="+",
        required=False,
        help="Fastq file, read1 of PE",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        nargs="+",
        required=False,
        help="Fastq file, read2 of PE",
    )
    parser.add_argument(
        "--mut",
        nargs="+",
        dest="mut",
        default=None,
        help="keyword of mutant fastq file, auto-find read1/2",
    )
    parser.add_argument(
        "--wt",
        nargs="+",
        dest="wt",
        default=None,
        help="keyword of wt fastq file, auto-find read1/2",
    )
    return parser


def get_args_index1(parser):
    # priority: index_list > extra_index > ...
    # parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
    #     help='Directory saving results, default: [pwd]')
    parser.add_argument(
        "-g",
        "--genome",
        default=None,
        help="The name of the genome, eg: dm6, hg38, default: [None]",
    )
    parser.add_argument(
        "--genome-index",
        dest="genome_index",
        default=None,
        help="align index of genome",
    )
    parser.add_argument(
        "-k",
        "--spikein",
        default=None,
        help="Spike-in genome, eg: dm6, hg38, default: [None]",
    )
    parser.add_argument(
        "--spikein-index",
        dest="spikein_index",
        default=None,
        help="align index of spikein",
    )
    parser.add_argument(
        "-x",
        "--extra-index",
        dest="extra_index",
        help="extra index, plasmid,te,piRNA_cluster, ... if specified, ignore -g, -k",
    )  # ???
    parser.add_argument(
        "--index-list",
        dest="index_list",
        nargs="+",
        help="A list of index, priority 1, ignore -x, -g, -k, ...",
    )
    parser.add_argument(
        "--to-rRNA", dest="to_rRNA", action="store_true", help="Align to rRNA"
    )
    parser.add_argument(
        "--rRNA-index",
        dest="rRNA_index",
        default=None,
        help="align index of rRNA",
    )
    # parser.add_argument('--to-chrM', action='store_true', dest='to_chrM',
    #     help='Align reads to mito-chromosome first')
    # parser.add_argument('--to-MT-trRNA', action='store_true', dest='to_MT_trRNA',
    #     help='Align reads to chrM, tRNA and rRNAs first')
    return parser


def get_args_opt(parser):
    # optional arguments - 1
    parser.add_argument(
        "--trimmed",
        action="store_true",
        help="Skip trimming, input reads are already trimmed",
    )
    parser.add_argument(
        "--cut",
        action="store_true",
        help="Cut reads to 50nt, equal to: --cut-to-length 50 --recursive",
    )
    parser.add_argument(
        "--cut-to-length",
        dest="cut_to_length",
        default=0,
        type=int,
        help="cut reads to specific length from tail, default: [0]",
    )
    parser.add_argument(
        "--recursive", action="store_true", help="trim adapter recursively"
    )
    # optional further
    parser.add_argument(
        "--binSize",
        dest="bin_size",
        default=50,
        type=int,
        help="binSize for downstream analysis, default [50]",
    )
    parser.add_argument(
        "--keep-dup",
        dest="rm_dup",
        action="store_false",
        help="keep duplicates",
    )
    parser.add_argument(
        "--keep-temp",
        dest="keep_temp",
        action="store_false",
        help="keep temp files",
    )
    parser.add_argument(
        "--gene-bed",
        dest="gene_bed",
        default=None,
        help="The BED or GTF of genes, for TSS enrichment analysis",
    )
    parser.add_argument(
        "--extend-reads",
        dest="extend_reads",
        action="store_true",
        help="Extend reads for bigwig",
    )
    parser.add_argument(
        "--center-reads",
        dest="center_reads",
        action="store_true",
        help="Center reads for bigwig",
    )
    parser.add_argument(
        "-p",
        "--threads",
        default=1,
        type=int,
        help="Number of threads to launch, default [1]",
    )
    parser.add_argument(
        "-j",
        "--parallel-jobs",
        dest="parallel_jobs",
        default=1,
        type=int,
        help="Number of jobs run in parallel, default: [1]",
    )
    parser.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        help="if specified, overwrite exists file",
    )
    parser.add_argument(
        "-f",
        "--fast",
        dest="fast_mode",
        action="store_true",
        help="Run pipeline in fast-mode",
    )
    return parser
