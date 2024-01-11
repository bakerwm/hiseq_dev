#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
The default files for ATAC
: x, is the out_dir
: smp_name, 
: fq1
: fq2
"""

import os


def get_atac_dirs(x, smp_name):
    """
    x is the out_dir
    """
    d = {
        "config_dir": "config",
        "raw_dir": "raw_data",
        "clean_dir": "clean_data",
        "align_dir": "align",
        "spikein_dir": "spikein",
        "bam_dir": "bam_files",
        "bg_dir": "bg_files",
        "bw_dir": "bw_files",
        "peak_dir": "peak",
        "motif_dir": "motif",
        "qc_dir": "qc",
        "report_dir": "report",
    }
    return {k: os.path.join(x, smp_name, v) for k, v in d.items()}


# def get_atac_files1(x, smp_name, fq1, fq2):
#     d = atac_basic_files(x, smp_name, fq1, fq2)
#     d.update(atac_qc_files(x, smp_name))
#     return d


def get_atac_files(x, smp_name, fq1, fq2):
    """
    x, is the out_dir
    """
    dd = get_atac_dirs(x, smp_name)  #
    config_dir = dd.get("config_dir", x)
    clean_dir = dd.get("clean_dir", x)
    bam_dir = dd.get("bam_dir", x)
    bw_dir = dd.get("bw_dir", x)
    peak_dir = dd.get("peak_dir", x)
    report_dir = dd.get("report_dir", x)
    d = {
        "config_yaml": os.path.join(config_dir, "config.yaml"),
        "trim_stat": os.path.join(clean_dir, smp_name + ".trim.stat"),
        "trim_json": os.path.join(clean_dir, smp_name + ".trim.json"),
        "align_scale_json": os.path.join(bam_dir, "scale.json"),
        "pcr_dup_json": os.path.join(bam_dir, "pcr_dup.json"),
        "bam_raw": os.path.join(bam_dir, smp_name + ".raw.bam"),
        "bam": os.path.join(bam_dir, smp_name + ".bam"),
        "bw": os.path.join(bw_dir, smp_name, smp_name + ".bigWig"),
        "peak": os.path.join(peak_dir, smp_name + "_peaks.narrowPeak"),
        "peak_seacr": os.path.join(peak_dir, smp_name + ".stringent.bed"),
        "peak_seacr_top001": os.path.join(
            peak_dir, smp_name + ".top0.01.stringent.bed"
        ),
        "report_log": os.path.join(report_dir, "report.log"),
        "report_html": os.path.join(report_dir, "HiSeq_report.html"),
    }
    d.update(atac_align_files(x, smp_name))
    d.update(atac_spikein_files(x, smp_name))
    if isinstance(fq1, str) and isinstance(fq2, str):
        d.update(atac_fq_files(x, smp_name, fq1, fq2))
    d.update(atac_qc_files(x, smp_name))
    return d


def atac_align_files(x, smp_name):
    """
    x, is the out_dir
    """
    dd = get_atac_dirs(x, smp_name)
    align_dir = dd.get("align_dir", x)
    d = {
        "align_bam": ".bam",
        "align_stat": ".align.stat",
        "align_json": ".align.json",
        "align_flagstat": ".align.flagstat",
        "unmap1": ".unmap.1.fastq",
        "unmap2": ".unmap.2.fastq",
    }
    return {k: os.path.join(align_dir, smp_name + v) for k, v in d.items()}


def atac_spikein_files(x, smp_name):
    """
    x is the spikein_dir
    """
    dd = get_atac_dirs(x, smp_name)
    spikein_dir = dd.get("spikein_dir", x)
    d = {
        "spikein_bam": ".bam",
        "spikein_scale_json": ".scale.json",
        "spikein_stat": ".align.stat",
        "spikein_json": ".align.json",
        "spikein_flagstat": ".align.flagstat",
        "spikein_unmap1": ".unmap.1.fastq",
        "spikein_unmap2": ".unmap.2.fastq",
    }
    return {k: os.path.join(spikein_dir, smp_name + v) for k, v in d.items()}


def atac_fq_files(x, smp_name, fq1, fq2):
    """
    x is the out_dir
    """
    dd = get_atac_dirs(x, smp_name)  #
    raw_dir = dd.get("raw_dir", x)
    clean_dir = dd.get("clean_dir", x)
    d = {
        "raw_fq1": os.path.join(raw_dir, os.path.basename(fq1)),
        "raw_fq2": os.path.join(raw_dir, os.path.basename(fq2)),
        "clean_fq1": os.path.join(clean_dir, os.path.basename(fq1)),
        "clean_fq2": os.path.join(clean_dir, os.path.basename(fq2)),
    }
    d.update(
        {
            "raw_fq_list": [d.get("raw_fq1"), d.get("raw_fq2")],
            "clean_fq_list": [d.get("clean_fq1"), d.get("clean_fq2")],
        }
    )
    return d


def atac_qc_files(x, smp_name):
    """
    x is the out_dir
    """
    dd = get_atac_dirs(x, smp_name)  #
    qc_dir = dd.get("qc_dir", x)
    d = {
        "trim_summary_json": "00.trim_summary.json",
        "align_summary_json": "01.alignment_summary.json",
        "dup_summary_json": "01.pcr_dup_summary.json",
        "lendist_csv": "02.length_distribution.fragsize.csv",
        "lendist_txt": "02.length_distribution.txt",
        "lendist_pdf": "02.length_distribution.fragsize.pdf",
        "frip_json": "03.FRiP.json",
        "tss_enrich_matrix": "04.tss_enrich.mat.gz",
        "tss_enrich_matrix_log": "04.tss_enrich.log",
        "tss_enrich_png": "04.tss_enrich.png",
        "tss_enrich_cmd": "04.tss_enrich.cmd.sh",
        "genebody_enrich_matrix": "05.genebody_enrich.mat.gz",
        "genebody_enrich_matrix_log": "05.genebody_enrich.log",
        "genebody_enrich_png": "05.genebody_enrich.png",
        "genebody_enrich_cmd": "05.genebody_enrich.cmd.sh",
        "bam_cor_npz": "06.bam_cor.npz",
        "bam_cor_counts": "06.bam_cor.counts.tab",
        "bam_cor_heatmap_png": "06.bam_cor.cor_heatmap.png",
        "bam_cor_pca_png": "06.bam_cor.cor_PCA.png",
        "peak_idr_png": "07.peak_idr.png",
        "peak_idr_txt": "07.peak_idr.txt",
        "peak_overlap_png": "08.peak_overlap.png",
        "peak_overlap_tiff": "08.peak_overlap.tiff",
        "bam_fingerprint_png": "09.fingerprint.png",
    }
    return {k: os.path.join(qc_dir, v) for k, v in d.items()}


# def main():
#     x = 'aaa'
#     a = atac_files('aaa', 'demo', 'pe_1.fq.gz', 'pe_2.fq.gz')
#     for k,v in a.items():
#         # print(v)
#         print('{}: {}'.format(k,v))


# if __name__ == '__main__':
#     main()
