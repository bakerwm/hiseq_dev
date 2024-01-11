#!/usr/bin/env python3
# -*- encoding:utf-8 -*-

"""
NAME: hiseq_p7.py 
DESCEIPTION:
what: Check the P7 sequence for Illumina Nextera and TruSeq libraries
how: 1. extract p7, at the 3' end of of read1 or read2; 
2. check regular i7-index (6-base or 8-base); 
3. also check the inline barcode from the beginning of read1, read2

# 1. Regular TruSeq and Nextera library
read1: --{insert}--{p7a}-{i7}-{p7b}
read2: --{insert}--{p5a}-{i5}-{p5b}

TruSeq:
p7a: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
i7:  6-base
p7b: ATCTCGTATGCCGTCTTCTGCTTG
p5a: AGATCGGAAGAGCGTCGTGTAGGGAAAGA
i5:  6-base
p5b: GTGTAGATCTCGGTGGTCGCCGTATCATT

Nextera:
p7a: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
i7:  8-base 
p7b: ATCTCGTATGCCGTCTTCTGCTTG
p5a: CTGTCTCTTATACACATCTGACGCTGCCGACGA
i5:  8-base
p5b: GTGTAGATCTCGGTGGTCGCCGTATCATT

# 2. Illumina library with inline barcodes
read1: -{bc/UMI}-{insert}--{p7a}-{i7}-{p7b}
read2: -{bc/UMI}-{insert}--{p5a}-{i5}-{p5b}

version-1: barcode 7-bp at the beginning of read2, eg: CLIP, 7-bp
version-2: UMI at the beginning of read1, eg: Phillip Zamore small RNA, 15-bp

================================================================================
Example output (json): top-10 for each group
{
    "name": "RNAseq_shWhite_0-2h_rep1_1",
    "count": 16313570,
    "p7a": "Truseq",
    "p5a": "Truseq",
    "i7": {
        "p7b_01": {
            "name": "p7b_01",
            "label": "TruSeq_p7b_01",
            "type": "TruSeq",
            "seq": "ATCTCGTATGCCGTCTTCTGCTTG",
            "count": 6442
        },
        "p7b_02": {
            "name": "p7b_02",
            "label": "Nextera_p7b_01",
            "type": "Nextera",
            "seq": "ATCTCGTATGCCGTCTTCTGCTTG",
            "count": 4031
        }
    },
    "i7": {
        "i7_01": {
            "name": "TruSeq_Index1",
            "label": "TruSeq_Index1",
            "type": "TruSeq",
            "seq": "ATCACG",
            "count": 8331
        },
        "i7_02": {
            "name": "Next_Ad2.1",
            "label": "Next_Ad2.1",
            "type": "Nextera",
            "seq": "TAAGGCGA",
            "count": 4401
        }
    },
    "p5a": {
        ""
    },
    "p5b": {
        ""
    },
    "i5": {
        ""
    }
}
================================================================================
history:
v4 - 2023-06-01
1. re-construct the script

v3 - 2021-06-15
1. guess library, 
2. support NSR

v-2 - 2021-05-11
1. support Nextera 
2. support large insert, (partial P7)

v-1 - 2021-05-11 
1. check full version of p7a, p7b
2. search barcode and i7


## to-do
1. merge guess_p5, guess_p7
2. reprot library structure
"""


import os
import sys
import re
import argparse
import pathlib
import shutil
import tempfile
import pyfastx
import logging
import hiseq  # for report
import Levenshtein as lev  # distance
from collections import Counter
from multiprocessing import Pool

# from hiseq.utils.fastx import Fastx
from hiseq.utils.utils import update_obj, log, run_shell_cmd, Config
from hiseq.utils.file import (
    check_dir,
    list_file,
    file_abspath,
    file_exists,
    file_abspath,
)
from hiseq.utils.seq import fx_name, list_fx
from hiseq.demx.sample_sheet import HiSeqIndex  # check index name/seq


logging.basicConfig(
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)
log.setLevel("INFO")


## tmp functions
def str_distance(x, y, partial=True):
    if isinstance(x, str) and isinstance(y, str):
        try:
            if partial:
                x = x[: len(y)]
                y = y[: len(x)]
            out = lev.distance(x, y)
        except:
            out = -1  # huge number
    else:
        out = -1
    return out


class HiseqAdapter(object):
    """
    support: TruSeq library (122nt+INS)
    # p5 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' # 58nt
    # p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' # 34nt
    # p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
    # p7 # 6nt

    support: Nextera library (128nt+INS)
    # p5 = 'AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' # 62nt
    # p7a = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # 34nt
    # p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
    # p7 # 8nt

    Sturcture
    P5a{i5}P5b--insert--P7a{i7}P7b

    TruSeq: P5a{58nt} + P7a{i7}P7b{64nt} = 122nt
    Nextera: P5a{62nt} + p7a{i7}P7b{66nt} = 128nt

    Current: single-index library (i7)
    """

    def __init__(self, lib="TruSeq"):
        lib = lib.lower()
        if not lib in ["truseq", "nextera"]:
            raise ValueError(
                "unknown lib, expect [truseq, nextera], got: {}".format(lib)
            )

    def get_adapters(self):
        #               read1->p7a->i7->p7b
        # p5b<-i5<-p5a<-read2
        ad = {
            "truseq": {
                "p7": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                "p7a": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                "p7b": "ATCTCGTATGCCGTCTTCTGCTTG",
                "p5": "AGATCGGAAGAGCGTCGTGTAGGGAAAGANNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT",
                "p5a": "AGATCGGAAGAGCGTCGTGTAGGGAAAGA",
                "p5b": "GTGTAGATCTCGGTGGTCGCCGTATCATT",
            },
            "nextera": {
                "p7": "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                "p7a": "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                "p7b": "ATCTCGTATGCCGTCTTCTGCTTG",
                "p5": "CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT",
                "p5a": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
                "p5b": "GTGTAGATCTCGGTGGTCGCCGTATCATT",
            },
        }
        return ad

    def get_adapter(self, lib="truseq", group="p7a"):
        ad = self.get_adapters()
        try:
            a = ad.get(lib, None)
            if isinstance(a, dict):
                out = a.get(group, None)
            else:
                out = None
        except:
            log.error(f"unknown lib={lib} and group={group}")
            out = None
        return out

    def guess_ad(self, x):
        ad = self.get_adatpers()
        out = None
        for lib in ad:  # truseq, nextera
            for k, v in ad[lib].items():
                if x == v:
                    out = f"{lib}:{k}"
        return out


def top_seq(x, n=3, rc=False, cutoff=0.05):
    out = {}
    if isinstance(x, list):
        t = len(x)
        a = Counter(x).most_common(n)
        for i in a:
            k, v = i
            rp = rev_comp(k) if rc else k
            rn = HiSeqIndex(rp).name  # NULL if not found
            if rn is None or rn == "NULL":  # in case, toml, 'NULL'
                rn = rp
            # seq, seq-rev, seq-name, pct, count
            if v / t > cutoff:
                out.update(
                    {
                        k: {
                            "seq": k,
                            "seq_revcomp": rp,
                            "name": rn,
                            "count": v,
                            "pct": f"{v/t:.4f}",
                        }
                    }
                )
        # sorted
        out = dict(
            sorted(out.items(), key=lambda x: x[1]["count"], reverse=True)
        )
    return out


def rev_comp(x):
    """Reverse complement DNAseq"""
    t = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
    }
    s = x[::-1]  # reverse
    return "".join([t.get(i, i) for i in s])  # complement


def guess_p7(x, cutoff=0.9, top_n=10000, out_json=None):
    ## bc -- p7a -- i7 -- p7b
    ad = HiseqAdapter().get_adapters()
    ## p7a
    p7at = HiseqAdapter().get_adapter("truseq", "p7a")
    p7an = HiseqAdapter().get_adapter("nextera", "p7a")
    pp1t = re.compile("([ACGTN]{6})" + p7at[:6])  # first 6-base
    pp1n = re.compile("([ACGTN]{6})" + p7an[:6])  # first 6-base
    ## p7a + p7b
    p7b = HiseqAdapter().get_adapter("truseq", "p7b")
    pp2t = re.compile(p7at[-6:] + "([ACGTN]{0,10})" + p7b[:6])  #
    pp2n = re.compile(p7an[-6:] + "([ACGTN]{0,10})" + p7b[:6])  #
    ## p7b (truseq == nextera)
    pp3 = re.compile(p7b[:6])  # first 6-base
    ## statistics
    i = 0  # total
    ## for NSR
    n_ct = 0
    n_ga = 0
    ## for barcode, first-7 base
    bc_list = []

    ## for TruSeq/Nextera
    n_p7a = 0
    n_p7b = 0
    n_truseq = 0
    n_nextera = 0
    is_nsr = 0  # 0=no, 1=yes
    i7_list = []
    for _, s, _, _ in pyfastx.Fastx(x):
        ## guess NSR BEGIN ##
        n_ct += int(s.startswith("CT"))
        n_ga += int(s.startswith("GA"))
        ## guess NSR END ##

        ## guess TruSeq/Nextera P7a ##
        s1t = pp1t.search(s)
        s1n = pp1n.search(s)
        # TruSeq
        if s1t:
            n_truseq += 1
            n_p7a += 1
            bc_list.append(s1t.group(1))  # barcode
            s2t = pp2t.search(s)
            if s2t:
                i7_list.append(s2t.group(1))  # i7 index
                n_p7b += 1
            else:
                pass
        # Nextera
        elif s1n:
            n_nextera += 1
            n_p7a += 1
            bc_list.append(s1n.group(1))  # barcode
            s2n = pp2n.search(s)
            if s2n:
                i7_list.append(s2n.group(1))  # i7 index
                n_p7b += 1
            else:
                pass
        else:
            pass
        # stop
        i += 1  # total count
        if i >= top_n:
            break

    ## summary
    ### 1. NSR
    if n_ct / i > cutoff:
        is_nsr = 1
        nsr = "read1"
    elif n_ga / i > cutoff:
        is_nsr = 1
        nsr = "read2"
    else:
        nsr = "Null"
    ### 2. P7
    top_i7 = top_seq(i7_list, n=3)
    is_truseq = int(n_truseq >= n_nextera and n_truseq > 0)
    is_nextera = int(n_nextera >= n_truseq and n_nextera > 0)
    top_bc = top_seq(bc_list, n=3, rc=True)  #

    # print('!A-1', n_truseq, n_nextera)
    out = {
        "name": os.path.basename(x),
        "is_truseq": is_truseq,
        "is_nextera": is_nextera,
        "is_nsr": is_nsr,
        "nsr": nsr,
        "i7": top_i7,
        "barcode": top_bc,
    }
    # save output
    if isinstance(out_json, str):
        try:
            Config().dump(out, out_json)
        except:
            log.error(f"failed writting to file: {out_json}")

    # return
    return out


def guess_p5(x, cutoff=0.9, top_n=10000, out_json=None):
    ## bc -- p5a -- i5 -- p5b
    ad = HiseqAdapter().get_adapters()
    ## p5a
    p5at = HiseqAdapter().get_adapter("truseq", "p5a")
    p5an = HiseqAdapter().get_adapter("nextera", "p5a")
    pp1t = re.compile("([ACGTN]{6})" + p5at[:6])  # first 6-base
    pp1n = re.compile("([ACGTN]{6})" + p5an[:6])  # first 6-base
    ## p5a + p5b
    p5b = HiseqAdapter().get_adapter("truseq", "p5b")
    pp2t = re.compile(p5at[-6:] + "([ACGTN]{0,10})" + p5b[:6])  #
    pp2n = re.compile(p5an[-6:] + "([ACGTN]{0,10})" + p5b[:6])  #
    ## p5b (truseq == nextera)
    pp3 = re.compile(p5b[:6])  # first 6-base
    ## statistics
    i = 0  # total
    ## for NSR
    n_ct = 0
    n_ga = 0
    ## for barcode, first-7 base
    bc_list = []

    ## for TruSeq/Nextera
    n_p5a = 0
    n_p5b = 0
    n_truseq = 0
    n_nextera = 0
    is_nsr = 0  # 0=no, 1=yes
    i5_list = []
    for _, s, _, _ in pyfastx.Fastx(x):
        ## guess NSR BEGIN ##
        n_ct += int(s.startswith("CT"))
        n_ga += int(s.startswith("GA"))
        ## guess NSR END ##

        ## guess TruSeq/Nextera p5a ##
        s1t = pp1t.search(s)
        s1n = pp1n.search(s)
        # TruSeq
        if s1t:
            n_truseq += 1
            n_p5a += 1
            bc_list.append(s1t.group(1))  # barcode
            s2t = pp2t.search(s)
            if s2t:
                i5_list.append(s2t.group(1))  # i5 index
                n_p5b += 1
            else:
                pass
        # Nextera
        elif s1n:
            n_nextera += 1
            n_p5a += 1
            bc_list.append(s1n.group(1))  # barcode
            s2n = pp2n.search(s)
            if s2n:
                i5_list.append(s2n.group(1))  # i5 index
                n_p5b += 1
            else:
                pass
        else:
            pass
        # stop
        i += 1  # total count
        if i >= top_n:
            break

    ## summary
    ### 1. NSR
    if n_ct / i > cutoff:
        is_nsr = 1
        nsr = "read1"
    elif n_ga / i > cutoff:
        is_nsr = 1
        nsr = "read2"
    else:
        nsr = "Null"
    ### 2. p5
    top_i5 = top_seq(i5_list, n=3)
    is_truseq = int(n_truseq >= n_nextera and n_truseq > 0)
    is_nextera = int(n_nextera >= n_truseq and n_nextera > 0)
    top_bc = top_seq(bc_list, n=3, rc=True)  #

    # print('!A-1', n_truseq, n_nextera)
    out = {
        "name": os.path.basename(x),
        "is_truseq": is_truseq,
        "is_nextera": is_nextera,
        "is_nsr": is_nsr,
        "nsr": nsr,
        "i5": top_i5,
        "barcode": top_bc,
    }
    # save output
    if isinstance(out_json, str):
        try:
            Config().dump(out, out_json)
        except:
            log.error(f"failed writting to file: {out_json}")

    # return
    return out


class HiseqP7(object):
    """
    Parse p7 from read1
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "fq1": None,
            "fq2": None,
            "out_dir": None,
            "parallel_jobs": 1,
            "overwrite": False,
            "top_n": 100000,
        }
        self = update_obj(self, args_init, force=False)
        self.fq = self.fq1  # update fq
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.fq_list = self.init_fq(self.fq)

    def init_fq(self, fq):
        """
        Make sure fq is read1 of PE reads
        """
        out = []
        if isinstance(fq, str):
            fq = os.path.abspath(fq)
            if os.path.isdir(fq):
                out = list_fx(fq)
                out = [i for i in out if "_1.f" in i]  # rep1
            elif os.path.isfile(fq):
                fname = fx_name(fq, fix_pe=False)
                if fname.endswith("_1") and os.path.exists(fq):
                    out = [fq]
            else:
                log.warning(
                    "illegal fq, str expect, got {}".format(type(fq).__name__)
                )
        elif isinstance(fq, list):
            out = [i for k in fq for i in self.init_fq(k)]
        else:
            log.warning(
                "illegal fq, str or list expect, got {}".format(
                    type(fq).__name__
                )
            )
        return out

    def run_single(self, x):
        # x_json = os.path.splitext(x)[0]+'.p7.json' #
        x_name = os.path.basename(os.path.splitext(x)[0])
        x_json = os.path.join(self.out_dir, x_name + ".p7.json")
        if file_exists(x_json) and not self.overwrite:
            log.info("file exists: {}".format(x_json))
        else:
            out = guess_p7(x, top_n=self.top_n, out_json=x_json)

    def run_multi(self, x):
        if self.parallel_jobs > 1 and isinstance(x, list):
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, x)
        else:
            for i in x:
                self.run_single(i)

    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, "bin", "hiseq_p7_report.R")
        report_html = os.path.join(self.out_dir, "HiSeq_P7_report.html")
        stdout = os.path.join(self.out_dir, "report.stdout")
        stderr = os.path.join(self.out_dir, "report.stderr")
        cmd_file = os.path.join(self.out_dir, "cmd.sh")
        cmd = " ".join(
            [
                shutil.which("Rscript"),
                qc_reportR,
                self.out_dir,
                "1> {}".format(stdout),
                "2> {}".format(stderr),
            ]
        )
        with open(cmd_file, "wt") as w:
            w.write(cmd + "\n")
        if os.path.exists(report_html) and self.overwrite is False:
            log.info("file exists, skip generating html.")
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error("failed, generating html file")

    def run(self):
        check_dir(self.out_dir)
        self.run_multi(self.fq_list)
        # self.report()


class HiseqP5(object):
    """
    Parse p5 from read2
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "fq1": None,
            "fq2": None,
            "out_dir": None,
            "parallel_jobs": 1,
            "overwrite": False,
            "top_n": 100000,
        }
        self = update_obj(self, args_init, force=False)
        self.fq = self.fq1  # update fq
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.fq_list = self.init_fq(self.fq)

    def init_fq(self, fq):
        """
        Make sure fq is read2 of PE reads
        """
        out = []
        if isinstance(fq, str):
            fq = os.path.abspath(fq)
            if os.path.isdir(fq):
                out = list_fx(fq)
                out = [i for i in out if "_2.f" in i]  # rep1
            elif os.path.isfile(fq):
                fname = fx_name(fq, fix_pe=False)
                if fname.endswith("_2") and os.path.exists(fq):
                    out = [fq]
            else:
                log.warning(
                    "illegal fq, str expect, got {}".format(type(fq).__name__)
                )
        elif isinstance(fq, list):
            out = [i for k in fq for i in self.init_fq(k)]
        else:
            log.warning(
                "illegal fq, str or list expect, got {}".format(
                    type(fq).__name__
                )
            )
        return out

    def run_single(self, x):
        x_name = os.path.basename(os.path.splitext(x)[0])
        x_json = os.path.join(self.out_dir, x_name + ".p5.json")
        if file_exists(x_json) and not self.overwrite:
            log.info("file exists: {}".format(x_json))
        else:
            out = guess_p5(x, top_n=self.top_n, out_json=x_json)

    def run_multi(self, x):
        if self.parallel_jobs > 1 and isinstance(x, list):
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, x)
        else:
            for i in x:
                self.run_single(i)

    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, "bin", "hiseq_p7_report.R")
        report_html = os.path.join(self.out_dir, "HiSeq_P5_report.html")
        stdout = os.path.join(self.out_dir, "report.stdout")
        stderr = os.path.join(self.out_dir, "report.stderr")
        cmd_file = os.path.join(self.out_dir, "cmd.sh")
        cmd = " ".join(
            [
                shutil.which("Rscript"),
                qc_reportR,
                self.out_dir,
                "1> {}".format(stdout),
                "2> {}".format(stderr),
            ]
        )
        with open(cmd_file, "wt") as w:
            w.write(cmd + "\n")
        if os.path.exists(report_html) and self.overwrite is False:
            log.info("file exists, skip generating html.")
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error("failed, generating html file")

    def run(self):
        check_dir(self.out_dir)
        self.run_multi(self.fq_list)
        # self.report()


class HiseqSmrna(object):
    def __init__(self, **kwargs):
        pass


def get_args():
    parser = argparse.ArgumentParser(description="hiseq i7")
    #     parser.add_argument('-i', '--fq', nargs='+', required=True,
    #         help='read1 of Paired end reads')
    parser.add_argument(
        "-1",
        "--fq1",
        nargs="+",
        required=True,
        help="read1 of Paired end reads",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        nargs="+",
        required=False,
        default=None,
        help="read2 of Paired end reads, optional, for smRNA only",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        required=True,
        help="The directory to save results.",
    )
    parser.add_argument(
        "-n",
        "--top-n",
        dest="top_n",
        type=int,
        default=100000,
        help="Top N reads to process, default: [100000]",
    )
    parser.add_argument(
        "-j",
        "--parallel-jobs",
        default=1,
        type=int,
        dest="parallel_jobs",
        help="Number of jobs run in parallel, default: [1]",
    )
    parser.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        help="overwrite the exists files",
    )
    parser.add_argument(
        "--smRNA",
        action="store_true",
        help="for small RNA analysis, UMI+barcode",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.8,
        help="cutoff for matching UMI/barcode for smRNA",
    )
    parser.add_argument(
        "--debug",
        dest="verbose",
        action="store_true",
        help="Show log message in details",
    )
    return parser


def main():
    args = vars(get_args().parse_args())
    if args["smRNA"]:
        HiseqSmrna(**args).run()
    else:
        # for read1
        args["fq"] = args["fq1"]  # only for fq1
        HiseqP7(**args).run()
        # for read2
        args["fq"] = args["fq2"]  # only for fq1
        HiseqP5(**args).run()


if __name__ == "__main__":
    main()


#
