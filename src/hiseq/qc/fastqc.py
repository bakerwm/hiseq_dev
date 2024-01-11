#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
FastQC for fastq files, using Falco/FastQC
1. Falco [https://github.com/smithlabcode/falco] # version 1.2.1, 2022-09-15
2. FastQC [https://github.com/s-andrews/FastQC]  # version 0.11.9, 
"""

import os
import sys
import re
import argparse
from shutil import which
from multiprocessing import Pool
from hiseq.qc.read_fastqc import ReadFastQC
from hiseq.utils.seq import fx_name, list_fx
from hiseq.utils.utils import (
    update_obj,
    log,
    run_shell_cmd,
    Config,
    gen_random_string,
)
from hiseq.utils.file import (
    file_exists,
    file_abspath,
    file_abspath,
    fix_out_dir,
    copy_file,
    remove_file,
    remove_dir,
    check_dir,
)
from hiseq.report.hiseq_report import HiSeqRpt


class FastQCR1(object):
    # for single fastq files
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "fq": None,  # str
            "out_dir": None,
            "threads": 1,
            "parallel_jobs": 1,
            "overwrite": False,
            "falco": None,
            "fastqc": None,
            "nogroup": False,
            "skip_html": False,
        }
        self = update_obj(self, args_init, force=False)
        if self.is_fq(self.fq):
            self.fq = file_abspath(self.fq)
            self.name = fx_name(self.fq)
            self.out_txt = os.path.join(
                self.out_dir, self.name + "_fastqc_data.txt"
            )
            self.out_json = os.path.join(
                self.out_dir, self.name + "_fastqc_data.json"
            )
            self.out_html = os.path.join(
                self.out_dir, self.name + "_fastqc_report.html"
            )
            self.out_dir = fix_out_dir(self.out_dir)

    def is_fq(self, x):
        """
        str and endswith fastq/fq(.gz)
        """
        if isinstance(x, str):
            p = re.compile("\.f(ast)?q(\.gz)?", flags=re.IGNORECASE)
            out = isinstance(p.search(x), re.Match)
        else:
            out = False
        return out

    def is_tool(self, name):
        return which(name) is not None

    def run_falco(self):
        """
        options:
        --outdir          Create all output files in the specified
        --nogroup         Disable grouping of bases for reads >50bp
        -skip-report       Skip generating HTML file
        -skip-summary    Skip short summary
        --quiet           Do not print more run info
        ## expect output
        out_dir/fastqc_data.txt
        out_dir/fastqc_reprot.html
        """
        # tmp files
        tmp_out_dir = os.path.join(self.out_dir, "tmp_" + gen_random_string(6))
        tmp_txt = os.path.join(tmp_out_dir, "fastqc_data.txt")
        tmp_html = os.path.join(tmp_out_dir, "fastqc_report.html")
        # run cmd
        cmd = " ".join(
            [
                self.falco,
                "{}".format("--nogroup" if self.nogroup else ""),
                "{}".format("--skip-report" if self.skip_html else ""),
                "-skip-summary --quiet",
                "--outdir {}".format(tmp_out_dir),
                self.fq,
            ]
        )
        # out txt
        if os.path.exists(self.out_json) and not self.overwrite:
            print("file exists: {}".format(self.out_json))
        else:
            # save command
            cmd_txt = os.path.join(self.out_dir, f"{self.name}.falco.sh")
            with open(cmd_txt, "wt") as w:
                w.write(cmd + "\n")
            try:
                run_shell_cmd(cmd)
                if os.path.exists(tmp_txt):
                    copy_file(tmp_txt, self.out_txt)
                    remove_file(tmp_txt, ask=False)
                if os.path.exists(tmp_html):
                    copy_file(tmp_html, self.out_html)
                    remove_file(tmp_html, ask=False)
                # clear tmp dir
                if file_exists(tmp_out_dir):
                    remove_dir(tmp_out_dir, ask=False)
                ReadFastQC(self.out_txt).save_as(self.out_json)
            except:
                print("Failed to run Falco")

    def run_fastqc(self):
        """
        options:
        --outdir Create all output files in the specified output directory.
        --nogroup
        --extract
        --quiet
        ## expect output
        {filename}_fastqc/fastqc_data.txt
        {filename}_fastqc.html
        """
        tmp_out_dir = os.path.join(self.out_dir, "tmp_" + gen_random_string(6))
        tmp_dir = os.path.join(tmp_out_dir, self.name + "_fastqc")
        tmp_txt = os.path.join(tmp_dir, "fastqc_data.txt")
        tmp_html = tmp_dir + ".html"
        cmd = " ".join(
            [
                self.fastqc,
                "{}".format("--nogroup" if self.nogroup else ""),
                "--extract --quiet",
                "--outdir {}".format(tmp_out_dir),
                self.fq,
            ]
        )
        # out txt
        if os.path.exists(self.out_json) and not self.overwrite:
            print("file exists: {}".format(self.out_json))
        else:
            cmd_txt = os.path.join(self.out_dir, f"{self.name}.fastqc.sh")
            with open(cmd_txt, "wt") as w:
                w.write(cmd + "\n")
            try:
                check_dir(tmp_dir)
                run_shell_cmd(cmd)
                if os.path.exists(tmp_txt):
                    copy_file(tmp_txt, self.out_txt)
                if os.path.exists(tmp_html):
                    copy_file(tmp_html, self.out_html)
                if file_exists(tmp_dir):
                    remove_dir(tmp_dir, ask=False, check_empty=False)
                if file_exists(tmp_out_dir):
                    remove_dir(tmp_out_dir, ask=False, check_empty=False)
                ReadFastQC(self.out_txt).save_as(self.out_json)
            except:
                print("Failed to run FastQC")

    def pick_tool(self):
        if not file_exists(self.falco) and self.is_tool("falco"):
            self.falco = which("falco")
        if not file_exists(self.fastqc) and self.is_tool("fastqc"):
            self.fastqc = which("fastqc")
        # choose which falco or fastqc
        if file_exists(self.falco) and self.cmd in ["falco", "auto"]:
            out = self.run_falco
        elif file_exists(self.fastqc) and self.cmd in ["fastqc", "auto"]:
            out = self.run_fastqc
        else:
            out = None
            # out = self.run_falco if file_exists(self.falco) else \
            #     self.run_fastqc if file_exists(self.fastqc) else None
            # if out is None:
            msg = "\n".join(
                "=" * 80,
                "Required tool not found",
                "Try to install either of falco or fastqc",
                "see: https://github.com/smithlabcode/falco",
                "or: https://github.com/s-andrews/FastQC",
                "in brief:",
                "$ conda install -c bioconda falco" "or",
                "$ conda install -c bioconda fastqc",
                "=" * 80,
            )
            print(msg)
        # out = self.run_fastqc
        return out

    def run(self):
        if self.is_fq(self.fq):
            qc = self.pick_tool()
            if qc:
                qc()
        else:
            log.warning("not a fastq file: {}".format(self.fq))


class FastQC(object):
    # for single fastq files
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "fq": None,  # str
            "out_dir": None,
            "threads": 1,
            "parallel_jobs": 1,
            "overwrite": False,
            "falco": None,
            "fastqc": None,
            "nogroup": False,
            "skip_html": False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "qc_rn"
        self.fq = self.list_fq_files(self.fq)
        self.out_dir = fix_out_dir(self.out_dir)
        self.init_files()
        Config().dump(self.__dict__.copy(), self.config_yaml)

    def init_files(self):
        self.project_dir = self.out_dir
        self.config_dir = os.path.join(self.project_dir, "config")
        self.report_dir = os.path.join(self.out_dir, "report")
        self.config_yaml = os.path.join(self.config_dir, "config.yaml")
        [
            check_dir(i, create_dirs=True)
            for i in [self.config_dir, self.report_dir]
        ]

    def list_fq_files(self, x):
        out = []
        if isinstance(x, str):
            if os.path.exists(x):
                if os.path.isdir(x):
                    out = list_fx(x)
                elif os.path.isfile(x):
                    out = [x]
                else:
                    log.warning("file or path expected, got {}".format(x))
            else:
                log.warning("path not exists. {}".format(x))
        elif isinstance(x, list):
            out = [j for i in x for j in self.list_fq_files(i)]
        else:
            log.warning("str, list expected, got {}".format(type(x).__name__))
        return out

    def run_single_fq(self, x):
        args = self.__dict__.copy()
        args.update({"fq": x})
        FastQCR1(**args).run()

    def run(self):
        if self.parallel_jobs > 1 and len(self.fq) > 1:
            ## Pool() run in parallel
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fq, self.fq)
        else:
            [self.run_single_fq(fq) for fq in self.fq]
        # Generate report
        HiSeqRpt(self.out_dir, overwrite=self.overwrite).run()


def get_args():
    parser = argparse.ArgumentParser(description="hiseq qc, fastqc")
    parser.add_argument(
        "-i",
        "--fq",
        nargs="+",
        required=True,
        help="reads in FASTQ files, or directory contains fastq files",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="The directory to save results.",
    )
    parser.add_argument(
        "-p",
        "--threads",
        type=int,
        default=1,
        help="Number of threads for each job, default: [1]",
    )
    parser.add_argument(
        "-j",
        "--parallel-jobs",
        type=int,
        default=1,
        dest="parallel_jobs",
        help="Number of jobs run in parallel, default: [1]",
    )
    parser.add_argument(
        "-O", "--overwrite", action="store_true", help="Overwrite exists file"
    )
    parser.add_argument(
        "-g",
        "--nogroup",
        action="store_true",
        help="Disable grouping of bases for reads >50bp",
    )
    parser.add_argument(
        "-H",
        "--skip-html",
        dest="skip_html",
        action="store_true",
        help="Skip generating HTML file",
    )
    parser.add_argument(
        "-c",
        "--cmd",
        default="falco",
        choices=["falco", "fastqc", "auto"],
        help="Choose the FastQC tools, [falco, fastqc], default: [auto]",
    )
    parser.add_argument(
        "--falco",
        default="falco",
        help="Specify the path of Falco, default: [falco]",
    )
    parser.add_argument(
        "--fastqc",
        default="fastqc",
        help="Specify the path of FastQC, default: [fastqc]",
    )
    return parser


def main():
    args = vars(get_args().parse_args())
    FastQC(**args).run()


if __name__ == "__main__":
    main()

#
