#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Align reads to index, using hisat2

Basic usage: 

1. SE
hisat2 -S -x index in.fq > out.sam 2> out.log

2. PE
hisat2 -S -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
"""

import os
import sys
import re
import pathlib
import shutil
import argparse

# from Levenshtein import distance
from hiseq.utils.utils import update_obj, Config, log, run_shell_cmd
from hiseq.utils.file import (
    check_file,
    check_dir,
    file_exists,
    file_abspath,
    symlink_file,
    remove_file,
)
from hiseq.utils.seq import (
    Fastx,
    check_fx,
    check_fx_args,
    check_fx_paired,
    fx_name,
)
from hiseq.align.align_index import AlignIndex
from hiseq.align.align_args import get_args_align


def parse_hisat2(x):
    """
    Wrapper hisat2 log
    for PE: only proper paired reads
    SE:
    100000 reads; of these:
      100000 (100.00%) were unpaired; of these:
        49955 (49.95%) aligned 0 times
        44883 (44.88%) aligned exactly 1 time
        5162 (5.16%) aligned >1 times
    50.05% overall alignment rate
    PE:
    100000 reads; of these:
      100000 (100.00%) were paired; of these:
        56761 (56.76%) aligned concordantly 0 times
        38882 (38.88%) aligned concordantly exactly 1 time
        4357 (4.36%) aligned concordantly >1 times
        ----
        56761 pairs aligned concordantly 0 times; of these:
          425 (0.75%) aligned discordantly 1 time
        ----
        56336 pairs aligned 0 times concordantly or discordantly; of these:
          112672 mates make up the pairs; of these:
            103130 (91.53%) aligned 0 times
            8593 (7.63%) aligned exactly 1 time
            949 (0.84%) aligned >1 times
    48.44% overall alignment rate
    unique, multiple, unmap, map, total
    """
    total = 0
    mapped = 0
    unmapped = 0
    multi = 0
    unique = 0
    out = 0
    is_paired = False
    warn_chunkmbs = False
    if check_file(x, check_empty=True):
        # processing
        with open(x) as r:
            for line in r:
                n = line.strip().split()[0]
                if not re.match("^[0-9]+$", n):
                    continue
                n = eval(n)
                # parsing
                if "were paired; of these:" in line:
                    total = n
                    is_paired = True
                elif "aligned concordantly 0 times" in line:
                    unmapped = n
                elif "aligned 0 times" in line:
                    unmapped = n
                elif "aligned concordantly exactly 1 time" in line:
                    unique = n
                elif "aligned concordantly >1 times" in line:
                    multi = n
                elif "reads; of these" in line and not is_paired:
                    total = n
                elif "aligned exactly 1 time" in line and not is_paired:
                    unique = n
                elif "aligned >1 times" in line and not is_paired:
                    multi = n
                else:
                    pass
    else:
        log.error("file not exists, {}".format(x))
    # msg
    if warn_chunkmbs:
        log.warning(
            "{}\nset --chunkmbs 128, to fix the errors".format(warn_chunkmbs)
        )
    return {
        "total": total,
        "map": unique + multi,
        "unique": unique,
        "multi": multi,
        "unmap": unmapped,
    }


class Hisat2Config(object):
    """Check args, prepare files for hisat2
    arguments
    output files
    parser ?!
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "aligner": "hisat2",
            "fq1": None,
            "fq2": None,
            "out_dir": None,
            "index": None,
            "index_name": None,
            "smp_name": None,
            "extra_para": None,
            "smp_name": None,
            "threads": 1,
            "overwrite": False,
            "n_map": 0,
            "unique_only": False,
            "keep_tmp": False,
            "keep_unmap": True,
            "large_insert": False,
            "default_hisat2": False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "hisat2_r1"
        self.init_fx()
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        # index name
        if not AlignIndex(self.index, self.aligner).is_valid():
            raise ValueError("index not valid, {}".format(self.index))
        if self.index_name is None:
            self.index_name = AlignIndex(self.index).index_name()
        # update files
        self.init_files()

    def init_fx(self):
        """
        Make sure, fx
        1. str
        2. fq1 exists
        3. fq2 None or exists
        """
        # if not check_fx(self.fq1):
        #     raise ValueError('--fq1, not exists, or empty')
        # if self.fq2 is not None and not check_fx(self.fq2):
        #     raise ValueError('--fq2, not exists, or empty')
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError("--fq1, --fq2 faild, not properly paired")
        # format
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.fx_format = Fastx(self.fq1).format  # fasta/q
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # sample
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=True)
        self.rep_list = os.path.join(self.out_dir, self.smp_name)

    def init_files(self):
        self.project_dir = os.path.join(
            self.out_dir, self.smp_name, self.index_name
        )
        self.config_dir = os.path.join(self.project_dir, "config")
        # output files
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_files = {
            #             'project_dir': self.project_dir,
            "config_yaml": os.path.join(self.config_dir, "config.yaml"),
            "cmd_shell": os.path.join(self.project_dir, "cmd.sh"),
            "bam": prefix + ".bam",
            "sam": prefix + ".sam",
            "unmap": prefix + ".unmap." + self.fx_format,
            "unmap1": prefix + ".unmap.1." + self.fx_format,  #
            "unmap2": prefix + ".unmap.2." + self.fx_format,  #
            "align_log": prefix + ".align.log",
            "align_stat": prefix + ".align.stat",
            "align_json": prefix + ".align.json",
            "align_flagstat": prefix + ".align.flagstat",
        }
        self = update_obj(self, default_files, force=True)
        check_dir([self.project_dir, self.config_dir], create_dirs=True)


class Hisat2(object):
    """
    Alignment, using hisat2
    Single index, SE/PE
    1. SE
    hisat2 -x index in.fq > out.sam 2> out.log
    2. PE
    hisat2 -x index -1 r1.fq -2 r2.fq > out.sam 2> out.log
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_local = Hisat2Config(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)  # update
        self.aligner = "hisat2"  # force changed
        self.get_cmd()
        Config().dump(self.__dict__, self.config_yaml)

    def get_cmd(self):
        """
        The command line
        ########################
        # hisat2 unique reads #
        ########################
        filt by tag: YT:Z:CP
        YT:Z: String representing alignment type
        CP: Concordant; DP: Discordant; UP: Unpaired Mate; UU: Unpaired.
        answer-1: https://www.biostars.org/p/19283/#19292
        answer-2: https://www.biostars.org/p/133094/#133127

        -F 2048: suppress Supplementary alignments

        or set: -q 10 # works for hisat2

        n_map: -k x
        extra_para: 'extra'
        # example:
        """
        args_fmt = "-f" if self.fx_format == "fasta" else "-q"
        args_nmap = "-k {}".format(self.n_map) if self.n_map > 0 else ""
        args_common = "--no-mixed --no-discordant"
        args_extra = self.extra_para if self.extra_para else ""
        args_large_ins = "-X 2000" if self.large_insert else ""
        if self.default_hisat2:
            args_nmap = ""
            args_common = ""
            args_extra = ""
            args_large_ins = ""
        if self.is_paired:
            args_io = " ".join(
                [
                    "--un-conc {}".format(self.unmap),
                    "-1 {} -2 {}".format(self.fq1, self.fq2),
                ]
            )
        else:
            args_io = "--un {} -U {}".format(self.unmap, self.fq1)
        # command-line
        cmd_main = " ".join(
            [
                "{}".format(shutil.which("hisat2")),
                "--mm -p {}".format(self.threads),
                args_fmt,
                args_common,
                args_nmap,
                args_extra,
                args_large_ins,
                "-x {}".format(self.index),
                args_io,
                "1> {} 2> {}".format(self.sam, self.align_log),
            ]
        )
        # for unique
        if self.is_paired:
            if self.unique_only:
                cmd_unique = " ".join(
                    [
                        "&& samtools view -Sub -F 0x4 -F 2048 -q 10 -f 2 {}".format(
                            self.sam
                        )
                    ]
                )
            else:
                cmd_unique = " ".join(
                    [
                        "&& samtools view -Sub -F 0x4 -F 2048 {}".format(
                            self.sam
                        ),
                    ]
                )
        else:
            if self.unique_only:
                cmd_unique = " ".join(
                    [
                        "&& samtools view -Sub -F 0x4 -F 2048 -q 10 {}".format(
                            self.sam
                        )
                    ]
                )
            else:
                cmd_unique = " ".join(
                    [
                        "&& samtools view -Sub -F 0x4 -F 2048 {}".format(
                            self.sam
                        )
                    ]
                )
        # add cmd
        self.cmd = " ".join(
            [
                cmd_main,
                cmd_unique,
                "| samtools sort -@ {} -o {} -".format(self.threads, self.bam),
                "&& samtools index {}".format(self.bam),
                "&& samtools flagstat {} > {}".format(
                    self.bam, self.align_flagstat
                ),
            ]
        )

    def run(self):
        if file_exists(self.bam) and not self.overwrite:
            log.info("Hisat2() skipped, file exists: {}".format(self.bam))
        else:
            # save cmd
            with open(self.cmd_shell, "wt") as w:
                w.write(self.cmd + "\n")
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error("Hisat2() failed, check {}".format(self.align_log))
            # output file
            if not check_file(self.bam, check_empty=True):
                log.error("Hisat2() failed, check {}".format(self.align_log))
        # rename unmap files, unmap_1.fastq -> unmap.1.fastq
        if not self.is_paired:
            self.unmap1, self.unmap2 = (self.unmap, None)
        # log
        df = parse_hisat2(self.align_log)
        df.update(
            {
                "name": self.smp_name,
                "index": self.index_name,
                "unique_only": self.unique_only,
            }
        )
        Config().dump(df, self.align_json)
        # remove temp files
        del_list = [self.sam]
        if not self.keep_unmap:
            del_list.extend([self.unmap1, self.unmap2, self.unmap])
        # if not self.keep_tmp:
        remove_file(del_list, ask=False)
        if not file_exists(self.bam):
            log.error("Hisat2() failed, no bam output: {}".format(self.bam))
        return (self.bam, self.unmap1, self.unmap2)


def get_args_io():
    example = "\n".join(
        [
            "Examples:",
            "$ python hisat2.py -1 f1.fq -x genome -o output",
            "# add extra para",
            '$ python hisat2.py -1 f1.fq -2 f2.fq -x genome -o output -X "-X 2000"',
            "# unique reads, update index_name",
            "$ python hisat2.py -1 f1.fq -x genome -o output -u -in 01.genome",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="run_hisat",
        description="run hisat2 program",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-1",
        "--fq1",
        required=True,
        help="Fasta/q file, read1 of PE, or SE read",
    )
    parser.add_argument(
        "-2",
        "--fq2",
        required=False,
        default=None,
        help="Fasta/q file, read2 of PE, or SE read, optional",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        default=None,
        help="Directory saving results, default: [cwd]",
    )
    parser.add_argument(
        "-x", "--index", required=True, help="The alignment index for hisat2"
    )
    parser.add_argument(
        "-in",
        "--index-name",
        default=None,
        dest="index_name",
        help="The name of the index",
    )
    parser.add_argument(
        "-n",
        "--smp-name",
        default=None,
        dest="smp_name",
        help="The name of the sample",
    )
    return parser


def get_args():
    return get_args_align(get_args_io())


def main():
    args = vars(get_args().parse_args())
    Hisat2(**args).run()


if __name__ == "__main__":
    main()

# EOF
