#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Quality control matrix for ATACseq analysis
1. raw_data: > 20M
2. clean_data: > 20M (clean_pct, too-short)
3. align: map 20M (chrM, map%)
4. peaks
5. FRiP
6. peak_overlap (multi: _rn)
7. frag_size
8. tss_enrich # optional
9. genebody_enrich # optional
10. bam_fingerprint # optional
11. bam_cor
"""

import sys
import os
import re
import glob
import shutil
from hiseq.fragsize.fragsize import BamFragSizeR1
from hiseq.utils.file import (
    check_file,
    check_dir,
    copy_file,
    copy_dir,
    file_exists,
    file_nrows,
)
from hiseq.utils.utils import (
    log,
    Config,
    run_shell_cmd,
    find_longest_common_str,
    hash_string,
)
from hiseq.utils.hiseq_utils import read_hiseq, list_hiseq_file, is_hiseq_dir

from hiseq.utils.bam import Bam, Bam2cor, Bam2fingerprint
from hiseq.utils.bed import PeakIDR, BedOverlap, PeakFRiP
from hiseq.atac.atac_utils import get_mito_count


################################################################################
def qc_trim_summary(x, hiseq_type="r1"):
    if not is_hiseq_dir(x, hiseq_type):
        log.error(
            "qc_trim_summary() skipped, not a {} dir: {}".format(hiseq_type, x)
        )
        return None
    a = read_hiseq(x, hiseq_type)  # for general usage
    # check: r1, rn, ...
    if a.is_hiseq_r1:
        # init
        d = {
            "name": a.smp_name,
            "input": 1,
            "output": 1,
            "out_pct": 100.0,
            "rm_pct": 0,
        }
        # update
        if file_exists(a.trim_json):
            d1 = Config().load(a.trim_json)
            d1["rm_pct"] = 100.0 - d1.get("out_pct", 0)
            d.update(d1)
        # update more
        d.update(
            {
                "out_pct": float("{:.2f}".format(d["out_pct"])),
                "rm_pct": float("{:.2f}".format(d["rm_pct"])),
            }
        )
        # Config().dump(d, a.trim_summary_json)
    elif a.is_hiseq_rn:
        r1_list = list_hiseq_file(x, "trim_summary_json", "r1")
        d = {}
        if isinstance(r1_list, list) and len(r1_list) > 1:
            for i in r1_list:
                di = Config().load(i)
                d1 = {
                    k: d.get(k, 0) + di.get(k, 0) for k in ["input", "output"]
                }
                d.update(d1)
            # update name, pct
            d.update(
                {
                    "name": a.smp_name,
                    "out_pct": round(
                        d.get("output", 0) / d.get("input", 1) * 100, 2
                    ),
                }
            )
            d["rm_pct"] = round(100.0 - d.get("out_pct", 0), 2)
        elif isinstance(r1_list, str) and file_exists(r1_list):
            d = Config().load(r1_list)
        else:
            log.error("qc_trim_summary() skipped")
        # Config().dump(d, a.trim_summary_json)
    else:
        log.warning("qc_trim_summary() skipped")
        # return None
    Config().dump(d, a.trim_summary_json)
    return d


def qc_align_summary(x, hiseq_type="auto"):
    """
    Organize the alignment:
    add pcr_dup: a.pcr_dup_json
    output:
    name, total, map, unique, multi, spikein, rRNA, unmap, dup, nodup
    """
    if not is_hiseq_dir(x, hiseq_type):
        log.error(
            "qc_align_summary() skipped, not a {} dir: {}".format(
                hiseq_type, x
            )
        )
        return None
    a = read_hiseq(x, hiseq_type)  # for general usage
    if a.is_hiseq_r1:
        # spikein, rRNA, genome
        n_total = 0
        if file_exists(a.spikein_json):
            sp = Config().load(a.spikein_json)
            n_sp = sp.get("map", 0) if isinstance(sp, dict) else 0
            n_total = sp.get("total", 0) if isinstance(sp, dict) else 0
        else:
            n_sp = 0
        # pcr_dup
        if file_exists(a.pcr_dup_json):
            sd = Config().load(a.pcr_dup_json)
            n_dup = sd.get("dup", 0)
            n_nodup = sd.get("nodup", 0)
        else:
            n_dup = n_nodup = 0
        # name, total, map, unique, multi, unmap
        if file_exists(a.align_json):
            df = Config().load(a.align_json)
            if n_total == 0:
                n_total = df.get("total", 1)
            df.update(
                {
                    "total": n_total,
                    "spikein": n_sp,
                    "chrM": get_mito_count(a.bam),
                    "dup": n_dup,
                    "nodup": n_nodup,
                }
            )
            Config().dump(df, a.align_summary_json)
    elif a.is_hiseq_rn:
        r1 = list_hiseq_file(x, "align_summary_json", "r1")
        df = {}
        for i in r1:
            di = Config().load(i)
            df = {
                k: df.get(k, 0) + di.get(k, 0)
                for k, v in di.items()
                if type(v) == int
            }  # merge values
            d_str = {
                k: v
                for k, v in di.items()
                if type(v) == bool or type(v) == str
            }
            df.update(d_str)  # unique_only, index, name
            df["name"] = list_hiseq_file(x, "smp_name", "auto")  # update name
        Config().dump(df, a.align_summary_json)
    else:
        log.warning("qc_align_summary() skipped, no align_json")
        return None


def qc_lendist(x, hiseq_type="r1"):
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_lendist() failed, not a hiseq dir: {}".format(x))
        return None
    if os.path.exists(a.lendist_csv) and a.overwrite is False:
        log.info("qc_lendist() skipped: file exists: {}".format(a.lendist_csv))
    else:
        b = BamFragSizeR1(
            bam=a.bam,
            out_dir=a.qc_dir,
            labels=None,
            csv_file=a.lendist_csv,
        )
        try:
            b.run()
        except:
            log.error("qc_lendist() failed, {}".format(a.lendist_csv))


def qc_frip(x, hiseq_type="r1"):
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_frip() failed, not a hiseq dir: {}".format(x))
        return None
    qc_frip_dir = os.path.join(a.qc_dir, "frip_files")
    if check_file(a.frip_json, check_empty=True):
        log.info("qc_frip() skipped, file exists: {}".format(a.frip_json))
    else:
        try:
            # dict: index, total, n, frip
            s = PeakFRiP(
                peak=a.peak,
                bam=a.bam,
                method="featureCounts",  # bedtools
                out_dir=qc_frip_dir,
            ).run()
            # number of peaks
            s.update(
                {
                    "name": a.project_name,
                    "n_peaks": file_nrows(a.peak),
                }
            )
            Config().dump(s, a.frip_json)
        except:
            log.error("PeakFRiP() failed, see: {}".format(a.frip_json))


################################################################################
## function for tss,genebody enrich: for general hiseq purpose
##
## computematrix, plotProfile from deeptools
##
def qc_tss_enrich_tool(
    x, hiseq_type="r1", bw_type="r1", subcmd="reference-point", **kwargs
):
    """
    1. get bw list
    2. get labels
    3. get title
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_tss_enrich_tool() failed, not a hiseq dir: {}".format(x))
        return None
    bed = getattr(a, "gene_bed", None)
    if not file_exists(bed):
        log.error(
            "qc_tss_enrich_tool() failed, bed not exists: {}".format(bed)
        )
        return None
    arg_bed = "-R {}".format(bed)
    # default values
    # -b 2000 -a 2000 --binSize 10
    upstream = kwargs.get("upstream", 1000)
    downstream = kwargs.get("downstream", 1000)
    regionbody = kwargs.get("regionbody", 2000)
    binsize = 10  # force
    # -b 2000 -a 2000 --binSize 50
    arg_body = "-b {} -a {} --binSize {}".format(upstream, downstream, binsize)
    # add genebody for scale-regions
    # subcmd: scale-regions, reference-point
    if subcmd == "scale-regions":
        arg_body += " -m {}".format(regionbody)
    # r1
    bw = ""
    arg_bw = ""
    arg_label = ""
    arg_title = ""
    per_group = "--perGroup"
    # r1: atac/rnaseq/cnr...
    if a.is_hiseq_r1:
        arg_title = "--plotTitle {}".format(a.smp_name)
        if a.hiseq_type.startswith("rnaseq_"):
            bw = " ".join([a.bw_fwd, a.bw_rev])
            arg_label = "--samplesLabel {} {}".format("fwd", "rev")
        else:
            bw = a.bw
            arg_label = "--samplesLabel {}".format(a.smp_name)
        arg_bw = "-S {}".format(bw)
    elif a.is_hiseq_rn:
        arg_title = "--plotTitle {}".format(a.smp_name)
        if a.hiseq_type.startswith("rnaseq"):
            b1 = list_hiseq_file(x, "bw_fwd", "r1")
            b2 = list_hiseq_file(x, "bw_rev", "r1")
            n1 = list_hiseq_file(x, "smp_name", "rn")
            n2 = list_hiseq_file(x, "smp_name", "rn")
            n1 = [i.replace(a.smp_name + "_", "") for i in n1]
            n2 = [i.replace(a.smp_name + "_", "") for i in n2]
            n_list = [i + k for i in n1 + n2 for k in ["_fwd", "_rev"]]
            bw_list = [a.bw_fwd, a.bw_rev] + b1 + b2
            n_list = ["merge_fwd", "merge_rev"] + n_list
            bw = " ".join(bw_list)
            arg_label = "--samplesLabel {}".format(" ".join(n_list))
            arg_title = "--plotTitle {}".format(a.smp_name)
        else:  # atac, cnr, chip, ...
            bw_list = list_hiseq_file(x, "bw", "r1")
            smp_name = list_hiseq_file(x, "smp_name", "r1")  # multi
            smp_name = [i.replace(a.smp_name + "_", "") for i in smp_name]
            # add merge
            bw_list.insert(0, a.bw)
            smp_name.insert(0, a.smp_name)
            # prepare args
            bw = " ".join(bw_list)
            arg_label = "--samplesLabel {}".format(" ".join(smp_name))
        arg_bw = "-S {}".format(bw)
    elif a.is_hiseq_rx:
        arg_title = "--plotTitle {}".format(a.smp_name)
        if a.hiseq_type in ["rnaseq_rx"]:
            # mut, wt
            bw_list = [a.mut_bw_fwd, a.mut_bw_rev, a.wt_bw_fwd, a.wt_bw_rev]
            arg_bw = " ".join(bw_list)
            # shorter name
            s = find_longest_common_str(a.mut_name, a.wt_name)
            s1 = a.mut_name.replace(s, "")
            s2 = a.wt_name.replace(s, "")
            ss = "{}.vs.{}".format(s1, s2)
            smp_name = [s1 + ".fwd", s1 + ".rev", s2 + "fwd", s2 + ".rev"]
            arg_label = "--samplesLabel {}".format(" ".join(smp_name))
        elif a.hiseq_type in ["chipseq_rx", "cnr_rx", "cnt_rx"]:
            # ip, input
            bw_list = [a.ip_bw, a.input_bw, a.bw]
            bw = " ".join(bw_list)
            # shorter name
            s = find_longest_common_str(a.ip_name, a.input_name)
            s1 = a.ip_name.replace(s, "")
            s2 = a.input_name.replace(s, "")
            ss = "{}.vs.{}".format(s1, s2)
            smp_name = [s1, s2, ss]
            arg_label = "--samplesLabel {}".format(" ".join(smp_name))
        else:
            pass
        arg_bw = "-S {}".format(bw)
    else:
        log.error(
            "qc_genebody_enrich() failed, unknown hiseq dir: {}".format(x)
        )
    if bw == "":
        return None
    arg_ma = " ".join([arg_label, arg_body, arg_bw, arg_bed])
    arg_plot = " ".join([arg_title, per_group])
    # return arguments
    return (arg_ma, arg_plot)


def qc_tss_enrich(x, hiseq_type="r1", bw_type="r1", **kwargs):
    """
    Calculate the Genebody enrichment

    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['r1']

    bw_type:  str
        The hiseq type of bigWig file, options: ['r1', 'rn', 'rx']
        default: ['r1']

    $ computeMatrix referencepoint -b -R gene.bed -S in.bw -o mat.gz
    $ plotProfile -m mat.gz -o tss.png
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_tss_enrich() failed, not a hiseq dir: {}".format(x))
        return None
    # arg_bw, arg_label, arg_title, per_group
    args = qc_tss_enrich_tool(
        x, hiseq_type, bw_type, subcmd="reference-point", **kwargs
    )
    if args is None:
        log.error("qc_tss_enrich() failed")
        return None
    cmd = " ".join(
        [
            "{}".format(shutil.which("computeMatrix")),
            "reference-point",
            "--referencePoint TSS",
            "--sortRegions descend --skipZeros",
            args[0],  # -R -S -b -a --binSize
            "-o {}".format(a.tss_enrich_matrix),
            "-p {}".format(a.threads),
            "2> {}".format(a.tss_enrich_matrix_log),
            "&& {}".format(shutil.which("plotProfile")),
            "-m {}".format(a.tss_enrich_matrix),
            "-o {}".format(a.tss_enrich_png),
            args[1],  # --plotTitle --perGroup
            "--dpi 300",
        ]
    )
    if file_exists(a.tss_enrich_png) and not a.overwrite:
        log.info("qc_tss() skipped, file exists: {}".format(a.tss_enrich_png))
    else:
        if not file_exists(getattr(a, "gene_bed", None)):
            log.error("qc_tss() skipped, gene_bed not found")
        else:
            with open(a.tss_enrich_cmd, "wt") as w:
                w.write(cmd + "\n")
            try:
                run_shell_cmd(cmd)
            except:
                log.error(
                    "failed to run computeMatrix, see: {}".format(
                        a.tss_enrich_matrix_log
                    )
                )


def qc_genebody_enrich(x, hiseq_type="r1", bw_type="r1", **kwargs):
    """
    Calculate the Genebody enrichment

    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['r1']

    bw_type:  str
        The hiseq type of bigWig file, options: ['r1', 'rn', 'rx']
        default: ['r1']

    $ computeMatrix scale-regions -R gene.bed -S f.bigWig -o mat.gz
    $ plotProfile -m mat.gz -o gene_body.png
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_genebody_enrich() failed, not a hiseq dir: {}".format(x))
        return None
    # arg_bw, arg_label, arg_title, per_group
    args = qc_tss_enrich_tool(
        x, hiseq_type, bw_type, subcmd="scale-regions", **kwargs
    )
    if args is None:
        log.error("qc_genebody_enrich() failed")
        return None
    cmd = " ".join(
        [
            "{}".format(shutil.which("computeMatrix")),
            "scale-regions",
            "-o {}".format(a.genebody_enrich_matrix),
            "--sortRegions descend --skipZeros",
            args[0],  # -R -S -b -a -m --binSize
            "--smartLabels",
            "-p {}".format(a.threads),
            "2> {}".format(a.genebody_enrich_matrix_log),
            "&& {}".format(shutil.which("plotProfile")),
            "-m {}".format(a.genebody_enrich_matrix),
            "-o {}".format(a.genebody_enrich_png),
            args[1],  # --plotTitle --perGroup
            "--dpi 300",
        ]
    )
    if file_exists(a.genebody_enrich_png) and not a.overwrite:
        log.info(
            "qc_genebody_enrich() skipped, file exists: {}".format(
                a.genebody_enrich_png
            )
        )
    else:
        if not file_exists(getattr(a, "gene_bed", None)):
            log.error("qc_tss() skipped, gene_bed not found")
        else:
            with open(a.genebody_enrich_cmd, "wt") as w:
                w.write(cmd + "\n")
            try:
                run_shell_cmd(cmd)
            except:
                log.error(
                    "qc_genebody_enrich() failed, see: {}".format(
                        a.genebody_enrich_matrix_log
                    )
                )


################################################################################


################################################################################
def qc_bam_cor(x, hiseq_type="rn", bam_type="r1"):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    bam_type:  str
        The hiseq type of bam file, options: ['r1', 'rn', 'rx']
        default: ['r1']

    Compute correlation (pearson) between replicates
    window = 500bp

    eg:
    multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
        --outRawCounts *counts.tab -b bam
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_bam_cor() failed, not a hiseq dir: {}".format(x))
        return None
    bam_list = list_hiseq_file(x, "bam", bam_type)
    args = {
        "bam_list": bam_list,
        "out_dir": a.qc_dir,
        "prefix": "06.bam_cor",
        "threads": a.threads,
        "overwrite": a.overwrite,
        "binsize": 500,  # a.binsize,
    }
    if file_exists(a.bam_cor_heatmap_png) and not a.overwrite:
        log.info("qc_bam_cor() skipped, file exists")
    else:
        if all(file_exists(bam_list)):
            Bam2cor(**args).run()
        else:
            log.error("qc_bam_cor() failed, bam files not exists")


def qc_peak_idr(x, hiseq_type="rn", peak_type="r1"):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    peak_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_peak_idr() failed, not a hiseq dir: {}".format(x))
        return None
    peak_list = list_hiseq_file(x, "peak", peak_type)
    args = {
        "peak_list": peak_list,
        "out_dir": a.qc_dir,
        "prefix": "07.peak_idr",
        "overwrite": a.overwrite,
    }
    if file_exists(a.peak_idr_png) and not a.overwrite:
        log.info("qc_peak_idr() skipped, file exists")
    else:
        if all(file_exists(peak_list)):
            PeakIDR(**args).run()
        else:
            log.error("qc_peak_idr() failed, peak files not exists")


def qc_peak_overlap(x, hiseq_type="rn", peak_type="r1"):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    peak_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_peak_overlap() failed, not a hiseq dir: {}".format(x))
        return None
    peak_list = list_hiseq_file(x, "peak", peak_type)
    args = {
        "peak_list": peak_list,
        "out_dir": a.qc_dir,
        "prefix": "08.peak_overlap",
        "overwrite": a.overwrite,
    }
    if file_exists(a.peak_overlap_png) and not a.overwrite:
        log.info("qc_peak_overwrite() skipped, file exists")
    else:
        if all(file_exists(peak_list)):
            # takes too-long for peaks > 20k
            n_peaks = [file_nrows(i) > 20000 for i in peak_list]
            if all([i > 20000 for i in n_peaks]):
                log.warning(
                    "qc_peak_overlap() skipped, too many peaks (>20000): [{}]".format(
                        "\,".join(list(map(str, n_peaks)))
                    )
                )
            else:
                BedOverlap(**args).run()
        else:
            log.error("qc_peak_overlap() failed, peak files not exists")


def qc_bam_fingerprint(x, hiseq_type="rn", bam_type="r1"):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    bam_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type)  # for general usage
    if not a.is_hiseq:
        log.error("qc_bam_fingerprint() failed, not a hiseq dir: {}".format(x))
        return None
    bam_list = list_hiseq_file(x, "bam", bam_type)
    if isinstance(bam_list, str):
        bam_list = [bam_list]
    args = {
        "bam_list": bam_list,
        "out_dir": a.qc_dir,
        "prefix": "09.fingerprint",
        "threads": a.threads,
    }
    if file_exists(a.bam_fingerprint_png) and not a.overwrite:
        log.info("qc_bam_fingerprint() skipped, file exists")
    else:
        if all(file_exists(bam_list)):
            Bam2fingerprint(**args).run()
        else:
            log.error("qc_bam_fingerprint() failed, peak files not exists")


def copy_hiseq_qc(x, hiseq_type="_merge"):
    """
    Copy the following data to dest dir, for hiseq_merge module
    The SHA-256 value of project_dir (first 7-character) was used to mark duplicate smp_names

    - report_html, data/
    - trim_summary_json
    - align_summary_json
    - dup_summary_json # in align_summary_json
    - frip_json
    - lendist_csv
    - tss_enrich_png
    - genebody_enrich_png
    - peak_overlap_png
    - bam_cor_heatmap_png
    - bam_cor_pca_png

    Parameters
    x  :  str
        path to the hiseq_merge dir
    keys : str
        key name of the file, in config
    dest : str
        path to the directory, saving files
    """
    log.info("copy hiseq qc files ...")
    key_list = [
        "config_yaml",
        "report_html",
        "trim_summary_json",
        "align_summary_json",
        "dup_summary_json",
        "frip_json",
        "lendist_csv",
        "tss_enrich_png",
        "genebody_enrich_png",
        "peak_overlap_png",
        "bam_cor_heatmap_png",
        "bam_cor_pca_png",
    ]
    a = read_hiseq(x, hiseq_type)  # for general usage
    hiseq_list = list_hiseq_file(x, "hiseq_list", "auto")
    # copy files to project_dir/data/
    for s in hiseq_list:
        if not is_hiseq_dir(s):
            continue  # skip
        b = read_hiseq(s)
        src_name = list_hiseq_file(s, "smp_name", "auto")
        src_out_dir = list_hiseq_file(s, "out_dir", "auto")
        src_name = getattr(b, "smp_name", None)  #
        src_out_dir = getattr(b, "out_dir", None)  #
        src_tag = hash_string(src_out_dir)[:7]  # first-7-character
        for key in key_list:
            src_file = list_hiseq_file(s, key, b.hiseq_type)
            if src_file is None:
                #                 log.info('hiseq_dir missing {} in {}'.format(key, s))
                continue  # skip
            # convert str to list
            if isinstance(src_file, str):
                src_file = [src_file]
            for src in src_file:  # multiple files
                if src is None:
                    continue
                if not os.path.exists(src):
                    continue
                dest_dir = os.path.join(a.data_dir, src_name + "." + src_tag)
                check_dir(dest_dir)
                copy_file(src, dest_dir)  # report.html
                # if report/data
                src_data_dir = os.path.join(os.path.dirname(src), "data")
                dest_data_dir = os.path.join(dest_dir, "data")
                check_dir(dest_data_dir)
                if os.path.exists(src_data_dir) and not os.path.exists(
                    dest_data_dir
                ):
                    copy_dir(src_data_dir, dest_data_dir)
