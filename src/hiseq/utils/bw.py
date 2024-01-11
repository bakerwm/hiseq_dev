#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions for bigWig files
bigWig
"""

import os
import shutil
from hiseq.utils.utils import log, run_shell_cmd
from hiseq.utils.file import check_dir


def bw_compare(bw1, bw2, bw_out, operation="log2", **kwargs):
    """
    Compare two bigWig files: ip over input
    example:
    bigwigCompare -b1 bw1 -b2 bw2 --operation
    {log2, ratio, subtract, add, mean, reciprocal_ratio,
    first, second}
    -o out.bw
    """
    bin_size = kwargs.get("bin_size", 50)
    threads = kwargs.get("threads", 1)
    overwrite = kwargs.get("overwrite", False)
    cmd = " ".join(
        [
            "{}".format(shutil.which("bigwigCompare")),
            "--bigwig1 {} --bigwig2 {}".format(bw1, bw2),
            "--operation {}".format(operation),
            "--skipZeroOverZero",
            "--skipNAs",
            "--binSize {}".format(bin_size),
            "-p {}".format(threads),
            "-o {}".format(bw_out),
        ]
    )
    # savd cmd
    bw_out_dir = os.path.dirname(bw_out)
    check_dir(bw_out_dir)
    cmd_txt = os.path.join(bw_out_dir, "cmd.txt")
    with open(cmd_txt, "wt") as w:
        w.write(cmd + "\n")
    # run
    if os.path.exists(bw_out) and not overwrite:
        log.info("bwCompare() skipped, file exists: {}".format(bw_out))
    else:
        run_shell_cmd(cmd)
