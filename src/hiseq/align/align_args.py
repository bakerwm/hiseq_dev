#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import argparse


def get_args_io1():
    example = "\n".join(
        [
            "Examples:",
            "1. auto-detect adapters",
            "$ python trim_r1.py -1 fq1 -2 fq2 -o out_dir -m 20",
            "2. specific 3-apdapter",
            "$ python trim_r1.py -1 fq1 -2 fq2 -o out_dir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="trim",
        description="trim adapters",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        required=True,
        help="reads in FASTQ files, support (*.gz), 1 file.",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        default=None,
        help="The second file of pair-end reads",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="Directory saving results, default: [cwd]",
    )
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="The name of the sample",
    )
    return parser


def get_args_index1(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-x", "--index", required=False, help="The alignment index for aligner"
    )
    group.add_argument(
        "-g",
        "--genome",
        default=None,
        help="The name of the genome, eg: dm6, mm10, hg38",
    )
    return parser


def get_args_io2():
    example = "\n".join(
        [
            "Examples:",
            "1. auto-detect adapters",
            "$ python trim_r1.py -1 fq1 -2 fq2 -o out_dir -m 20",
            "2. specific 3-apdapter",
            "$ python trim_r1.py -1 fq1 -2 fq2 -o out_dir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="trim",
        description="trim adapters",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1",
        "--fq1",
        dest="fq1",
        nargs="+",
        required=True,
        help="reads in FASTQ files, support (*.gz), 1 file.",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        dest="fq2",
        nargs="+",
        default=None,
        help="The second file of pair-end reads",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="Directory saving results, default: [cwd]",
    )
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="The name of the sample",
    )
    return parser


def get_args_index2(parser):
    # priority: index_list > extra_index > ...
    parser.add_argument(
        "--index-list", dest="index_list", nargs="+", help="A list of index"
    )
    parser.add_argument(
        "-x",
        "--extra-index",
        dest="extra_index",
        help="Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k",
    )
    parser.add_argument(
        "-g",
        "--genome",
        default=None,
        help="The name of the genome, [dm6, hg38, mm10]",
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
        help="Spike-in genome : dm3, hg19, hg38, mm10, default: None",
    )
    parser.add_argument(
        "--spikein-index",
        dest="spikein_index",
        default=None,
        help="align index of spikein",
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
    parser.add_argument(
        "--to-chrM",
        action="store_true",
        dest="to_chrM",
        help="Align reads to mito-chromosome first",
    )
    parser.add_argument(
        "--to-MT-trRNA",
        action="store_true",
        dest="to_MT_trRNA",
        help="Align reads to chrM, tRNA and rRNAs first",
    )
    return parser


def get_args_align(parser):
    # extra para
    parser.add_argument(
        "-u",
        "--unique-only",
        action="store_true",
        dest="unique_only",
        help="Report unique mapped reads only",
    )
    parser.add_argument(
        "--n-map",
        dest="n_map",
        type=int,
        default=0,
        help="Number of hits per read, default [0]",
    )
    parser.add_argument(
        "-l",
        "--large-insert",
        action="store_true",
        dest="large_insert",
        help="For large insert, use: -X 1000 --chunkmbs 128",
    )
    parser.add_argument(
        "--keep-tmp",
        dest="keep_tmp",
        action="store_true",
        help="save temp files",
    )
    parser.add_argument(
        "--rm-unmap",
        dest="keep_unmap",
        action="store_false",
        help="remove unmap fastq files",
    )
    parser.add_argument(
        "-X",
        "--extra-para",
        dest="extra_para",
        default=None,
        help='Add extra parameters, eg: "-X 2000"',
    )
    parser.add_argument(
        "-p",
        "--threads",
        default=4,
        type=int,
        help="Number of threads, default: [4]",
    )
    parser.add_argument(
        "-j",
        "--parallel-jobs",
        default=1,
        type=int,
        help="Number of jobs to run in parallel, default: [1]",
    )
    parser.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        help="Overwrite the exist files",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Show message in details"
    )
    parser.add_argument(
        "-a",
        "--aligner",
        default="bowtie2",
        type=str,
        help="The aligner for alignment, default: [bowtie2]",
    )
    parser.add_argument(
        "--genome-dir",
        dest="genome_dir",
        type=str,
        default=None,
        help="The directory of genome files, default: [${HOME}/data/genome/]",
    )
    return parser
