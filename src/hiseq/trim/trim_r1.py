#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: single fx
1. Trim adapter,
2. remove pcr dup
3. trim n-bases from read
3.1 TruSeq (NSR)
    - cut 7,-7; from both ends of read (subseq)
3.2 TruSeq (iCLIP)
    - cut 9; from read1
3.3 TruSeq (eCLIP)
    - cut 10,-7 from read1
    - cut 7,-10 from read2
"""

import os
import re
import shutil
import pathlib
import argparse
import pandas as pd
import pyfastx
from xopen import xopen
from hiseq.utils.utils import log, update_obj, Config, get_date
from hiseq.trim.cutadapt import Cutadapt, get_args_trim
from hiseq.utils.seq import Fastx, check_fx, check_fx_paired, fx_name, readfq
from hiseq.utils.file import (
    check_dir,
    file_exists,
    file_abspath,
    copy_file,
    symlink_file,
    remove_file,
)
from hiseq.report.hiseq_report import HiSeqRpt


class TrimR1(object):
    """
    for single fastq file
    1. trim adapter
    2. remove dup
    3. cut n-bases
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        object = TrimR1Config(**self.__dict__.copy())
        self = update_obj(self, object.__dict__.copy(), force=True)
        Config().dump(self.__dict__.copy(), self.config_yaml)

    def get_msg(self):
        msg = "\n".join(
            [
                "-" * 80,
                "{:>20s} : {}".format("program", "Trim_r1"),
                "{:>20s} : {}".format("Date", get_date()),
                "{:>20s} : {}".format("config", self.config_yaml),
                "{:>20s} : {}".format("fq1", self.fq1),
                "{:>20s} : {}".format("fq2", self.fq2),
                "{:>20s} : {}".format("out_dir", self.out_dir),
                "{:>20s} : {}".format("library_type", self.library_type),
                "{:>20s} : {}".format("len_min", self.len_min),
                "{:>20s} : {}".format("len_max", self.len_max),
                "{:>20s} : {}".format("rm_dup", self.rm_dup),
                "{:>20s} : {}".format("cut_before_trim", self.cut_before_trim),
                "{:>20s} : {}".format("cut_after_trim", self.cut_after_trim),
                "{:>20s} : {}".format("3-adapter (read1)", self.adapter3),
                "{:>20s} : {}".format("3-adapter (read2)", self.Adapter3),
                "{:>20s} : {}".format("5-adapter (read1)", self.adapter5),
                "{:>20s} : {}".format("5-adapter (read2)", self.Adapter5),
                "{:>20s} : {}".format("remove poly(N)", self.rm_polyN),
                "{:>20s} : {}".format(
                    "remove untrimmed", self.discard_untrimmed
                ),
                "-" * 80,
            ]
        )
        return msg

    def run(self):
        print(self.get_msg())
        Trim = TrimR1pe if self.is_paired else TrimR1se
        Trim(**self.__dict__.copy()).run()
        HiSeqRpt(x=self.project_dir, overwrite=self.overwrite).run()


class TrimR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "library_type": None,  # auto
            "fq1": None,  # str
            "fq2": None,  # str
            "out_dir": None,  # str
            "smp_name": None,
            "threads": 4,
            "rm_dup": False,
            "cut_before_trim": "0",
            "cut_after_trim": "0",  # skip, str
            "keep_tmp": False,
            "len_min": 0,
            "len_max": 0,
            "adapter3": None,
            "Adapter3": None,
            "adapter5": None,
            "Adapter5": None,
            "rm_polyN": False,
            "discard_untrimmed": False,
            "overwrite": False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "trim_r1"  #
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        # self.init_fq()
        # fix cut_after_trim
        if self.cut_after_trim is None:
            self.cut_after_trim = "0"  # default
        self.cut_after_trim = str(self.cut_after_trim)
        self.cut2_l, self.cut2_r = self.parse_cut_arg(self.cut_after_trim)
        self.cut2_skipped = self.cut2_l + abs(self.cut2_r) == 0

    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.out_dir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, "config")
        self.log_dir = os.path.join(self.project_dir, "log")
        self.tmp_dir = os.path.join(self.project_dir, "tmp")
        self.cutadapt_dir = os.path.join(self.tmp_dir, "01_cutadapt")
        self.rm_dup_dir = os.path.join(self.tmp_dir, "02_rm_dup")
        self.cut2_dir = os.path.join(self.tmp_dir, "03_cut_after_trim")
        # files
        default_files = {
            "config_yaml": os.path.join(self.config_dir, "config.yaml"),
            "rm_dup_fq": os.path.join(
                self.rm_dup_dir, self.smp_name + ".fq.gz"
            ),
            "rm_dup_fq1": os.path.join(
                self.rm_dup_dir, self.smp_name + "_1.fq.gz"
            ),
            "rm_dup_fq2": os.path.join(
                self.rm_dup_dir, self.smp_name + "_2.fq.gz"
            ),
            "cut2_fq": os.path.join(self.cut2_dir, self.smp_name + ".fq.gz"),
            "cut2_fq1": os.path.join(
                self.cut2_dir, self.smp_name + "_1.fq.gz"
            ),
            "cut2_fq2": os.path.join(
                self.cut2_dir, self.smp_name + "_2.fq.gz"
            ),
            "clean_fq": os.path.join(
                self.project_dir, self.smp_name + ".fq.gz"
            ),
            "clean_fq1": os.path.join(
                self.project_dir, self.smp_name + "_1.fq.gz"
            ),
            "clean_fq2": os.path.join(
                self.project_dir, self.smp_name + "_2.fq.gz"
            ),
            "trim_stat": os.path.join(
                self.project_dir, self.smp_name + ".trim.stat"
            ),
            "trim_json": os.path.join(
                self.project_dir, self.smp_name + ".trim.json"
            ),
            "trim_yaml": os.path.join(
                self.project_dir, self.smp_name + ".trim.yaml"
            ),
        }
        self = update_obj(self, default_files, force=True)  # key
        check_dir(
            [
                self.config_dir,
                self.cutadapt_dir,
                self.rm_dup_dir,
                self.cut2_dir,
                self.log_dir,
            ]
        )

    def is_valid_cut(self, s):
        if isinstance(s, str):
            s = re.sub("[^\d,-]", "", s)  # sanitize string
            p = re.compile("(^\d+)(,(-\d+))?$", re.IGNORECASE)
            out = isinstance(p.match(s), re.Match)
        else:
            out = False
        return out

    def parse_cut_arg(self, s):
        """
        for 1-indexed
        cut2: '10', '9,-9', '5', ...
        '10': [10:]
        '0,-9': [:-9]
        '9, -9': [9:-9]
        make sure: valid
        return: start, end
        """
        if self.is_valid_cut(s):
            p = re.compile("(^\d+)(,(-\d+))?$", re.IGNORECASE)
            m = p.match(s)
            l = eval(m.group(1))
            r = eval(m.group(3)) if m.group(2) else 0
        else:
            raise ValueError(
                "unknown cut arg: {}, {}".format(s, type(s).__name__)
            )
        return [l, r]


class TrimR1se(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_local = self.__dict__.copy()
        local_obj = TrimR1Config(**args_local)
        self = update_obj(self, local_obj, force=True)
        # Config().dump(self.__dict__.copy(), self.config_yaml)

    def cutadapt(self):
        args = self.__dict__.copy()
        args.update({"out_dir": self.cutadapt_dir})
        # cut1_args = self.__dict__.copy()
        # cut1_args['out_dir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**args)
        cut1.run()  #
        cut1_fq = cut1.clean_fq
        n_total, n_cut1, _ = cut1.parse_log()
        copy_file(cut1.log, self.log_dir)  # save log
        return (cut1_fq, n_total, n_cut1)

    def rm_pcr_dup(self, fq_in, fq_out):
        if check_fx(fq_in) and isinstance(fq_out, str):
            if self.rm_dup:
                if file_exists(fq_out) and not self.overwrite:
                    log.info("rm_dup() skipped, file exists")
                elif file_exists(fq_in):
                    Fastx(fq_in).collapse(fq_out, fq_out=True)
                else:
                    log.error("dup() failed: {}".format(fq_in))
                n_rm_dup = Fastx(self.rm_dup_fq).number_of_seq()
            else:
                symlink_file(fq_in, fq_out, absolute_path=False)
                n_rm_dup = -1
            return n_rm_dup

    def cut2(self, fq_in, fq_out):
        """
        sub-function
        make sure
        s: str
        l: int, >= 0
        r: int, <= 0
        """

        def substr(s, l, r):
            s1 = s[l:]  # trim left
            if r < 0:
                s1 = s1[:r]  # trim right
            return s1

        # run
        n_cut2_rm = 0
        if self.cut2_skipped:
            symlink_file(fq_in, fq_out, absolute_path=False)
        else:
            try:
                with xopen(fq_out, "wt") as w, xopen(fq_in) as r:
                    min_len = sum(
                        [self.len_min, self.cut2_l, abs(self.cut2_r)]
                    )
                    for name, seq, qual, cmt in readfq(r):
                        if len(seq) < min_len:
                            n_cut2_rm += 1
                            continue  # drop short seq
                        # update name
                        if cmt:
                            name = "{} {}".format(name, cmt)
                        # sub
                        seq = substr(seq, self.cut2_l, self.cut2_r)
                        qual = substr(qual, self.cut2_l, self.cut2_r)
                        # output
                        w.write("@{}\n{}\n+\n{}\n".format(name, seq, qual))
            except IOError as e:
                log.error(e)
        return n_cut2_rm

    def wrap(self, n_total, n_cut1, n_rm_dup, n_cut2):
        # 4. wrap files
        if n_total < 1:
            n_total = 1
        n_cut1_rm = n_total - n_cut1
        n_rm_dup_rm = n_cut1 - n_rm_dup
        n_cut2_rm = n_rm_dup - n_cut2
        n_clean_pct = "{:.2f}".format(n_cut2 / n_total * 100)
        # header
        header = [
            "#name",
            "total",
            "too_short",
            "dup",
            "too_short2",
            "clean",
            "percent",
        ]
        s = [
            self.smp_name,
            n_total,
            n_cut1_rm,
            n_rm_dup_rm,
            n_cut2_rm,
            n_cut2,
            n_clean_pct,
        ]
        s = list(map(str, s))
        msg = "\t".join(header) + "\n" + "\t".join(s)
        with open(self.trim_stat, "wt") as w:
            w.write(msg + "\n")
        # save to yaml
        df_stat = {
            "name": self.smp_name,
            "total": int(n_total),
            "too_short": int(n_cut1_rm),
            "dup": int(n_rm_dup_rm),
            "too_short2": int(n_cut2_rm),
            "clean": int(n_cut2),
            "percent": float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # save fq files, cutadapt log
        if file_exists(self.cut2_fq) and not file_exists(self.clean_fq):
            copy_file(self.cut2_fq, self.clean_fq)

    def run(self):
        ## 1.cutadapt; 2.rm_dup; 3.cut2
        if file_exists(self.clean_fq) and not self.overwrite:
            log.info("Trim() skipped, file exists: {}".format(self.clean_fq))
            return None
        # 1. cutadapt
        cut1_fq, n_total, n_cut1 = self.cutadapt()
        # 2. remove PCR dup
        n_rm_dup = self.rm_pcr_dup(cut1_fq, self.rm_dup_fq)
        n_rm_dup = n_cut1 if n_rm_dup < 0 else n_rm_dup
        # 3. cut2, further adapters
        n_cut2_rm = self.cut2(self.rm_dup_fq, self.cut2_fq)
        n_cut2 = n_rm_dup - n_cut2_rm
        # 4. wrap files
        self.wrap(n_total, n_cut1, n_rm_dup, n_cut2)
        # 5. remove temp files
        del_list = [cut1_fq, self.rm_dup_fq, self.cut2_fq]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)


class TrimR1pe(object):
    """
    for single fastq file: Paired-End (PE)
    1. trim adapter
    2. remove dup
    3. cut n-bases
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        obj = TrimR1Config(**self.__dict__.copy())
        self = update_obj(self, obj, force=True)

    def cutadapt(self):
        args = self.__dict__.copy()
        args.update({"out_dir": self.cutadapt_dir})
        # cut1_args = self.__dict__.copy()
        # cut1_args['out_dir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**args)
        cut1.run()  #
        cut1_fq1, cut1_fq2 = [cut1.clean_fq1, cut1.clean_fq2]
        n_total, n_cut1, _ = cut1.parse_log()
        copy_file(cut1.log, self.log_dir)  # save log
        return (cut1_fq1, cut1_fq2, n_total, n_cut1)

    def rm_pcr_dup(self, fq1_in, fq2_in, fq1_out, fq2_out):
        log.info("rm_dup() skipped, only works in single-end mode")
        symlink_file(fq1_in, fq1_out, absolute_path=False)
        symlink_file(fq2_in, fq2_out, absolute_path=False)
        # n_rm_dup = -1 # n_cut1

    def cut2(self, fq1_in, fq2_in, fq1_out, fq2_out):
        """
        sub-function
        make sure
        s: str
        l: int, >= 0
        r: int, <= 0
        """

        def substr(s, l, r):
            s1 = s[l:]  # trim left
            if r < 0:
                s1 = s1[:r]  # trim right
            return s1

        n_cut2_rm = 0
        if self.cut2_skipped:
            log.info("cut_after_trim(), skipped")
            symlink_file(fq1_in, fq1_out, absolute_path=False)
            symlink_file(fq2_in, fq2_out, absolute_path=False)
        else:
            seq_min = sum([self.len_min, self.cut2_l, abs(self.cut2_r)])
            try:
                with xopen(fq1_out, "wt") as w1, xopen(
                    fq2_out, "wt"
                ) as w2, xopen(fq1_in) as r1, xopen(fq2_in) as r2:
                    for f1, f2 in zip(readfq(r1), readfq(r2)):
                        # name,seq,qual,comment
                        if len(f1[1]) < seq_min or len(f2[1]) < seq_min:
                            n_cut2_rm += 1
                            continue  # drop short seq
                        # update name
                        n1 = "{} {}".format(f1[0], f1[3]) if f1[3] else f1[0]
                        n2 = "{} {}".format(f2[0], f2[3]) if f2[3] else f2[0]
                        # sub
                        s1 = substr(f1[1], self.cut2_l, self.cut2_r)
                        q1 = substr(f1[2], self.cut2_l, self.cut2_r)
                        s2 = substr(f2[1], self.cut2_l, self.cut2_r)
                        q2 = substr(f2[2], self.cut2_l, self.cut2_r)
                        # output
                        w1.write("@{}\n{}\n+\n{}\n".format(n1, s1, q1))
                        w2.write("@{}\n{}\n+\n{}\n".format(n2, s2, q2))
            except IOError as e:
                log.error(e)
        return n_cut2_rm

    def wrap(self, n_total, n_cut1, n_rm_dup, n_cut2):
        if n_total < 1:
            n_total = 1
        n_cut1_rm = n_total - n_cut1
        n_rm_dup_rm = n_cut1 - n_rm_dup
        n_cut2_rm = n_rm_dup - n_cut2
        n_clean_pct = "{:.2f}".format(n_cut2 / n_total * 100)
        # header
        header = [
            "#name",
            "total",
            "too_short",
            "dup",
            "too_short2",
            "clean",
            "percent",
        ]
        s = [
            self.smp_name,
            n_total,
            n_cut1_rm,
            n_rm_dup_rm,
            n_cut2_rm,
            n_cut2,
            n_clean_pct,
        ]
        s = list(map(str, s))
        msg = "\t".join(header) + "\n" + "\t".join(s)
        with open(self.trim_stat, "wt") as w:
            w.write(msg + "\n")
        # save to yaml
        df_stat = {
            "name": self.smp_name,
            "total": int(n_total),
            "too_short": int(n_cut1_rm),
            "dup": int(n_rm_dup_rm),
            "too_short2": int(n_cut2_rm),
            "clean": int(n_cut2),
            "percent": float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        copy_file(self.cut2_fq1, self.clean_fq1)
        copy_file(self.cut2_fq2, self.clean_fq2)

    def run(self):
        ## re-run
        if (
            all(file_exists([self.clean_fq1, self.clean_fq2]))
            and not self.overwrite
        ):
            log.info(
                "Trim() skipped, file exists: {}, {}".format(
                    self.clean_fq1, self.clean_fq2
                )
            )
            return None
        # 1. cut adapter
        cut1_fq1, cut1_fq2, n_total, n_cut1 = self.cutadapt()
        # 2. remove PCR dup
        self.rm_pcr_dup(cut1_fq1, cut1_fq2, self.rm_dup_fq1, self.rm_dup_fq2)
        n_rm_dup = n_cut1  # rm_dup skipped for PE reads
        # 3. cut, further adapters
        n_cut2_rm = self.cut2(
            self.rm_dup_fq1, self.rm_dup_fq2, self.cut2_fq1, self.cut2_fq2
        )
        n_cut2 = n_rm_dup - n_cut2_rm
        # 4. wrap files
        self.wrap(n_total, n_cut1, n_rm_dup, n_cut2)
        # 5. remove temp files
        del_list = [
            cut1_fq1,
            cut1_fq2,
            self.rm_dup_fq1,
            self.rm_dup_fq2,
            self.cut2_fq1,
            self.cut2_fq2,
        ]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)


def get_args_i():
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
    return parser


def get_args_o(parser):
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="The directory to save results.",
    )
    parser.add_argument(
        "-n",
        "--smp-name",
        dest="smp_name",
        default=None,
        help="The prefix of output files",
    )
    parser.add_argument(
        "--rm-dup",
        dest="rm_dup",
        action="store_true",
        help="remove duplicates",
    )
    parser.add_argument(
        "--cut2",
        "--cut-after-trim",
        dest="cut_after_trim",
        default="0",
        help="cut after trimming; example: 10: cut 10 from left; 0,-6; cut 6 from right; 9,-9, cut 9 from both ends, default [0]",
    )
    parser.add_argument(
        "--save-tmp",
        dest="keep_tmp",
        action="store_true",
        help="Save temp files",
    )
    return parser


def get_args():
    parser = get_args_o(get_args_i())  # I/O
    parser = get_args_trim(parser)
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    # group.add_argument("-q", "--quiet", action="store_true")
    return parser


def main():
    args = vars(get_args().parse_args())
    TrimR1(**args).run()


if __name__ == "__main__":
    main()

#
