#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare bam files
correlation: PCA, cor, using deeptools
"""

# import os
# import pathlib
import argparse

# from shutil import which
# from hiseq.utils.file import file_exists, check_path, file_prefix
from hiseq.utils.utils import (
    log,
)  # , Config, run_shell_cmd, get_date, update_obj
from hiseq.utils.bam import Bam2cor  # , Bw2cor


def get_args():
    parser = argparse.ArgumentParser(
        description="hiseq bam2cor -i bam -o out_dir"
    )
    parser.add_argument(
        "-i",
        "--bam-list",
        dest="bam_list",
        nargs="+",
        required=True,
        help="BAM files or bigWig files",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="output directory to save results",
    )
    parser.add_argument(
        "-m",
        "--cor-method",
        default="pearson",
        choices=["pearson", "spearman"],
        help="method to calculate correlation, default: [pearson]",
    )
    parser.add_argument(
        "-np",
        "--no-plot",
        dest="no_plot",
        action="store_false",
        help="do not make plots",
    )
    parser.add_argument(
        "-n",
        "--prefix",
        default=None,
        help="set the prefix for output files, default: [multibam]",
    )
    parser.add_argument(
        "-p",
        "--threads",
        default=1,
        type=int,
        help="number of threads, default: [1]",
    )
    parser.add_argument(
        "-b",
        "--binsize",
        default=500,
        type=int,
        help="set binSize for bigWig, default: [500]",
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        action="store_true",
        help="Whether overwrite exists files",
    )
    return parser


def bam2cor(**kwargs):
    bam_list = kwargs.get("bam_list", [])
    if all([i.endswith(".bam") for i in bam_list]):
        print("Bam2cor")
        Bam2cor(**kwargs).run()
    else:
        log.error("no bam files found")


def main():
    args = vars(get_args().parse_args())
    bam2cor(**args)


if __name__ == "__main__":
    main()

#
