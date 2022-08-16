#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
General modules for ATACseq analysis

analysis-module:
"""

import os
# import sys
import re
# import glob
# import shutil
import pysam
from hiseq.trim.trim_r1 import TrimR1
from hiseq.align.align import Align
from hiseq.bam2bw.bam2bw import Bam2bw
from hiseq.callpeak.callpeak import CallPeak
from hiseq.utils.bam import Bam
from hiseq.utils.file import check_file, symlink_file, file_exists, list_dir
from hiseq.utils.utils import log, Config, run_shell_cmd
from hiseq.utils.hiseq_utils import read_hiseq, list_hiseq_file, is_hiseq_dir
from hiseq.align.align_index import check_index_args


"""
Main:
1. trim
2. align spikein, extra, ...
3. align genome
4. call peak
5. bam to bw (optional)
6. qc
"""
def hiseq_trim(x, hiseq_type='_r1'):
    """
    Parameters
    ----------
    x: str
        The path to AtacR1() dir, hiseq_type=atac_r1
    option-1: cut reads from 3' end, to X-nt in length (default: 50)
    operation:
    1. do trimming
    2. create symlink
    path -> (AtacReader) -> args
    fq1, fq2
    clean_fq1, clean_fq2, clean_dir
    cut_to_length, recursive
    """
    if not is_hiseq_dir(x, '_r1'):
        log.error('hiseq_atac_trim() skipped, not a {} dir: {}'.format(hiseq_type, x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # do-the-trimming
    fq1, fq2 = a.raw_fq_list
    clean_fq1, clean_fq2 = a.clean_fq_list
    # whether to trim or not
    args = {
        'fq1': fq1,
        'fq2': fq2,
        'out_dir': a.clean_dir,
        'smp_name': a.smp_name,
        'library_type': None,
        'len_min': 20,
        'cut_to_length': a.cut_to_length,
        'recursive': a.recursive,
        'parallel_jobs': 1 # do not allowed > 1 !!!!
    }
    if a.trimmed:
        symlink_file(fq1, clean_fq1)
        symlink_file(fq2, clean_fq2)
    else:
        trim = TrimR1(**args)
        trim.run()
        symlink_file(trim.clean_fq1, clean_fq1)
        symlink_file(trim.clean_fq2, clean_fq2)
        symlink_file(trim.trim_json, a.trim_json)


def hiseq_align_spikein(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to AtacR1() dir, hiseq_type=atac_r1
    Align reads to reference genome, using bowtie2
    --sensitive --local -X 2000  # --devotail
    samtools view -bhS -f 2 -F 1804
    exclude: -F
        4    0x4    Read unmapped,
        8    0x8    Mate unmapped,
        256  0x100  Not primary alignment,
        512  0x200  Reads fails platform quality checks,
        1024 0x400  Read is PCR or optical duplicate
    : do not rm_dup :
    """
    if not is_hiseq_dir(x, '_r1'):
        log.error('hiseq_align_spikein() skipped, not a hiseq_r1 dir: {}'.format(x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # init index: spikein, spikein_index
    args = {
        'aligner': a.aligner,
        'genome': None,
        'genome_index': None,
        'spikein': a.spikein,
        'spikein_index': a.spikein_index,
    }
    index_list = check_index_args(**args)
    # update for alignment
    if is_hiseq_dir(x, 'atac_'):
        align_extra = '-I 10 -X 2000'
    elif is_hiseq_dir(x, 'cn'):
        align_extra = '-I 10 -X 700'
    else:
        align_extra = '' # empty
    # update optional
    args_align = {
        'aligner': a.aligner,
        'fq1': a.clean_fq1,
        'fq2': a.clean_fq2,
        'out_dir': a.spikein_dir,
        'unique_only': False,
        'index_list': index_list,
        'genome': None,
        'genome_index': None,
        'spikein': None,
        'spikein_index': None,
        'to_rRNA': None,
        'rRNA_index': None,
        'extra_index': None,
        'max_fragment': 2000,
        'extra_para': align_extra,
    }
    if file_exists(a.spikein_bam) and not a.overwrite:
        log.info('hiseq_align_spikein() skipped, file exists: {}'.format(
            a.spikein_bam))
    else:
        Align(**args_align).run()
    # copy files; go to align_r1 directory
    d_list = list_dir(a.spikein_dir, include_dir=True)
    d_list = [i for i in d_list if os.path.isdir(i)]
    for i in d_list:
        t = read_hiseq(i, 'alignment_rn') # alignment_rn
        if t.is_hiseq:
            symlink_file(t.bam, a.spikein_bam)
            symlink_file(t.align_json, a.spikein_json)
            symlink_file(t.align_flagstat, a.spikein_flagstat)
            symlink_file(t.unmap1, a.unmap1)
            symlink_file(t.unmap2, a.unmap2)
    # calculate norm scale
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=True)


def hiseq_align_genome(x, hiseq_type='_r1'):
    """
    Parameters
    ----------
    x: str
        The path to HiseqR1() dir, hiseq_type=hiseq_r1
    Align reads to reference genome, using bowtie2
    --sensitive --local -I 10 -X 700  # --no-mixed --no-discordant --devotail
    samtools view -bhS -f 2 -F 1804
    ######
    cnr: "-I 10 -X 700"
    atac: "-I 10 -X 2000"
    exclude: -F
        4    0x4    Read unmapped,
        8    0x8    Mate unmapped,
        256  0x100  Not primary alignment,
        512  0x200  Reads fails platform quality checks,
        1024 0x400  Read is PCR or optical duplicate
    """
    if not is_hiseq_dir(x, '_r1'):
        log.error('hiseq_align_genome() skipped, not a hiseq_r1 dir: {}'.format(x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # init index: spikein, spikein_index
    args = {
        'aligner': a.aligner,
        'genome': a.genome,
        'genome_index': a.genome_index,
        'spikein': None,
        'spikein_index': None,
        'to_MT': False, # on genome
        'to_chrM': False, # on genome
        'to_te': True, # false, if available !!! force
        'to_piRC': True, # false, if available !!! force
        'extra_index': a.extra_index,
        'index_list': a.index_list, # optional ???
    }
    index_list = check_index_args(**args)
    # update for alignment
    if is_hiseq_dir(x, 'atac_'):
        align_extra = '-I 10 -X 2000'
    elif is_hiseq_dir(x, 'cn'):
        align_extra = '-I 10 -X 700'
    else:
        align_extra = '' # empty
    # update optional
    # !!!; fq1,fq2 from unmap of spikein?
    args_align = {
        'aligner': a.aligner,
        'fq1': a.clean_fq1,
        'fq2': a.clean_fq2,
        'out_dir': a.align_dir,
        'unique_only': True,
        'index_list': index_list,
        'genome': a.genome, # !!!
        'genome_index': None,
        'spikein': None,
        'spikein_index': None,
        'to_rRNA': None,
        'rRNA_index': None,
        'extra_index': None,
        'max_fragment': 2000,
        'extra_para': align_extra,
    }
    if file_exists(a.spikein_bam) and not a.overwrite:
        log.info('hiseq_align_spikein() skipped, file exists: {}'.format(
            a.spikein_bam))
    else:
        Align(**args_align).run()
    # copy files; go to align_r1 directory
    d = os.path.join(a.align_dir, a.smp_name)
    if is_hiseq_dir(d, '_rn'): # fq=1, index=N
        dd = read_hiseq(d)
        symlink_file(dd.bam, a.bam_raw)
        symlink_file(dd.align_stat, a.align_stat)
        symlink_file(dd.align_json, a.align_json)
        symlink_file(dd.align_flagstat, a.align_flagstat)
    # remove PCR dup
    if a.rm_dup:
        if file_exists(a.bam_raw):
            if check_file(a.bam, check_empty=True) and not a.overwrite:
                log.info('rm_dup() skipped, file exists: {}'.format(
                    a.bam))
            else:
                Bam(a.bam_raw).rm_dup(outfile=a.bam)
    else:
        symlink_file(a.bam_raw, a.bam)
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=False)


def hiseq_merge_trim(x, hiseq_type='_rn'):
    """
    only works for atac_rn
    """
    if not is_hiseq_dir(x, '_rn'):
        log.error('hiseq_merge_bam() skipped, not a hiseq_rn dir: {}'.format(x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # ['clean', 'dup', 'name', 'percent', 'too_short', 'too_short2', 'total']
    tx = list_hiseq_file(x, 'trim_json', 'r1')
    d = {}
    if isinstance(tx, str):
        tx = [tx]
    for i in tx:
        di = Config().load(i)
        for k,v in di.items():
            vv = v if isinstance(v, str) else d.get(k, 0) + v
            d.update({k:vv})
    d.update({
        'name': a.smp_name,
        'percent': round(d.get('clean', 0) / d.get('total', 1) * 100, 1),
    })
    Config().dump(d, a.trim_json)


def hiseq_merge_bam(x, hiseq_type='_rn'):
    """
    only works for atac_rn
    """
    if not is_hiseq_dir(x, '_rn'):
        log.error('hiseq_merge_bam() skipped, not a hiseq_rn dir: {}'.format(x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    bam_list = list_hiseq_file(x, 'bam', 'r1') # single / multiple bam files
    if isinstance(bam_list, list) and len(bam_list) > 1:
        # 1. multiple bam files
        if not all(file_exists(bam_list)):
            log.error('could not merge bam files: {}'.format(a.bam_raw))
            return None
        # run
        if file_exists(a.bam_raw) and not a.overwrite:
            log.info('merge_bam() skipped, file exists: {}'.format(a.bam_raw))
        else:
            cmd = ' '.join([
                'samtools merge -',
                ' '.join(bam_list),
                '| samtools sort -o {} -'.format(a.bam_raw),
                '&& samtools index {}'.format(a.bam_raw)])
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')
        # rm_dup
        if a.rm_dup:
            if file_exists(a.bam_raw):
                if check_file(a.bam, check_empty=True) and not a.overwrite:
                    log.info('rm_dup() skipped, file exists: {}'.format(
                        a.bam))
                else:
                    Bam(a.bam_raw).rm_dup(a.bam)
        else:
            symlink_file(a.bam_raw, a.bam)
        # check-point
        if not file_exists(a.bam):
            log.error('could not merge bam files: {}'.format(a.bam_raw))
            return None
        # update align_stat/json; missing: (align_stat, align_flagstat)
        d = {'name': a.smp_name}
        for aj in list_hiseq_file(x, 'align_json', 'r1'):
            da = Config().load(aj)
            d.update({
                'index': da.get('index', None),
                'unique_only': da.get('unique_only', False),
                'total': d.get('total', 0) + da.get('total', 0),
                'map': d.get('map', 0) + da.get('map', 0),
                'unique': d.get('unique', 0) + da.get('unique', 0),
                'multi': d.get('multi', 0) + da.get('multi', 0)
            })
        Config().dump(d, a.align_json)
        # calculate norm scale
        s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=False) # to genome
    elif isinstance(bam_list, str):
        # 2. single bam file
        log.info('merge() skipped, Only 1 replicate detected')
        dd = read_hiseq(bam_list)
        tx = ['bam_raw', 'bam', 'align_stat', 'align_json', 'align_flagstat', 'align_scale_json']
        for i in tx:
            src = list_hiseq_file(bam_list, i, '_r1')
            dest = getattr(a, i)
            symlink_file(src, dest)


def hiseq_copy_r1(x, hiseq_type='_rn'):
    """
    copy files from rep1
    """
    if not is_hiseq_dir(x, '_rn'):
        log.error('hiseq_merge_bam() skipped, not a hiseq_rn dir: {}'.format(x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    bam_list = list_hiseq_file(x, 'bam', 'r1') # single / multiple bam files
    if isinstance(bam_list, str):
        log.warning('merge() skipped, Only 1 replicate detected')
        k_list = [
            'bam', 'bw', 'peak', 'peak_seacr', 'peak_seacr_top001',
            'align_scale_json', 'trim_json', 'align_json', 'pcr_dup_json'
        ]
        # get the rep list
        rep_list =  list_hiseq_file(x, 'rep_list', '_rn')
        r1_dir = rep_list[0]
        # rep_dir = self.rep_list[0] # first one
        for k in k_list:
            k_from = list_hiseq_file(r1_dir, k, 'r1')
            k_to = list_hiseq_file(x, k, 'rn')
            symlink_file(k_from, k_to) # update list_hiseq_file
        # copy bam.bai
        r1_bam_bai = list_hiseq_file(r1_dir, 'bam', 'r1')+'.bai'
        m_bam_bai = list_hiseq_file(x, 'bam', 'rn')+'.bai'
        symlink_file(r1_bam_bai, m_bam_bai)
        # copy qc
        m_qc_dir = list_hiseq_file(x, 'qc_dir', 'rn') # dest
        r1_qc_dir = list_hiseq_file(r1_dir, 'qc_dir', 'r1')
        r1_qc_files = list_dir(r1_qc_dir, include_dir=True) # update list_hiseq_file
        for f in r1_qc_files:
            symlink_file(f, m_qc_dir) # to qc_dir
        # update: bam index
        Bam(self.bam).index()


def hiseq_norm_scale(x, hiseq_type='_r1', by_spikein=False, norm=1000000):
    """
    Parameters
    ----------
    x: str
        The path to CnrRn() dir, hiseq_type=cnr_rn
    cal the norm scale: combine rep_list
    (fragments in bam files: reads, paired-reads)
    bam.count
    """
    if not is_hiseq_dir(x, hiseq_type):
        log.error('hiseq_norm_scale() skipped, not a {} dir: {}'.format(hiseq_type, x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # choose norm
    by_spikein = by_spikein and file_exists(a.spikein_bam)
    bam = a.spikein_bam if by_spikein  else a.bam
    if file_exists(a.align_scale_json):
        out = Config().load(a.align_scale_json)
    else:
        n = Bam(bam).count(reads=False) # paired
        try:
            s = round(norm / n, 4)
        except ZeroDivisionError as e:
            log.error(e)
            s = 1.0
        # output
        out = {
            'smp_name': a.smp_name,
            'is_spikein': by_spikein,
            'norm': norm,
            'map': n,
            'scale': s
        }
        Config().dump(out, a.align_scale_json)
    return out


def hiseq_pcr_dup(x, hiseq_type='r1'):
    """
    parse PCR-dup: r1, rn
    dup = align_json (bam_raw) - align_scale_json (bam)
    """
    if not is_hiseq_dir(x, hiseq_type):
        log.error('hiseq_merge_bam() skipped, not a {} dir: {}'.format(hiseq_type, x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    raw = list_hiseq_file(x, 'align_json', hiseq_type)
    nodup = list_hiseq_file(x, 'align_scale_json', hiseq_type)
    if all(file_exists([raw, nodup])):
        d1 = Config().load(raw)
        d2 = Config().load(nodup)
        d0 = {
            'name': d1.get('name'),
            'total': d1.get('total', 0),
            'nodup': d2.get('map', 0),
            'dup': d1.get('map', 0) - d2.get('map', 0),
        }
        d0['dup_pct'] = round(d0['dup'] / d0['total'], 4) if d0['total'] > 0 else 0
        Config().dump(d0, a.pcr_dup_json)
        return d0


def get_mito_count(x):
    """
    Count reads on chrM (MT), ...
    on genome, search for chrM, chrMT, ...
    """
    try:
        lines = pysam.idxstats(x).splitlines()
        lines = [i for i in lines if re.search('^(chrM|MT)', i)]
        n_mt = sum(
            map(int, [x.split("\t")[2]
                      for x in lines if not x.startswith("#")]))
    except IndexError as e:
        msg = "cannot read bam file, msg=%s, data=%s".format(
            (e, lines))
        log.error(msg)
        n_mt = 0
    return n_mt


def hiseq_call_peak(x, hiseq_type='_r1'):
    """
    Call peaks using MACS2 and SEACR
    -f BAMPE
    """
    if not is_hiseq_dir(x, hiseq_type):
        log.error('hiseq_call_peak() skipped, not a {} dir: {}'.format(hiseq_type, x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    # choose ip,input bam files
    if is_hiseq_dir(x, '_rx'):
        ip_bam = list_hiseq_file(x, 'ip_bam', hiseq_type)
        input_bam = list_hiseq_file(x, 'input_bam', hiseq_type)
    else:
        ip_bam = list_hiseq_file(x, 'bam', hiseq_type)
        input_bam = None
    # cmd
    args = {
        'ip': ip_bam, # rm_dup
        'input': input_bam,
        'out_dir': a.peak_dir,
        'prefix': a.smp_name,
        'genome': a.genome,
        'genome_size': None, # a.genome_size,
        'genome_size_file': a.genome_size_file
    }
    CallPeak(method='macs2', **args).run()
    CallPeak(method='seacr', **args).run()


def hiseq_bam2bw(x, hiseq_type='_r1'):
    """
    hiseq_bam_to_bw()
    Check the scale only for : TE/Genome
    normalizeUsing: RPGC, CPM
    """
    if not is_hiseq_dir(x, hiseq_type):
        log.error('hiseq_bam2bw() skipped, not a {} dir: {}'.format(hiseq_type, x))
        return None
    a = read_hiseq(x, hiseq_type) # for general usage
    strand_specific = a.hiseq_type.startswith('rna')
    args = {
        'bam': a.bam,
        'prefix': a.smp_name,
        'out_dir': a.bw_dir,
        'binSize': a.bin_size, # default: 50
        'strand_specific': strand_specific, # rna_ only
        'genome': a.genome,
        'scaleFactor': 1.0, # force
        'normalizeUsing': 'RPGC',
        'overwrite': a.overwrite,
        'genome_size': None, # a.genome_size,
        'extendReads': True, # a.extend_reads, # force !!!
        'centerReads': True, # a.center_reads, # force !!!
    }
    Bam2bw(**args).run()


def hiseq_bw_compare(x, hiseq_type='_rx'):
    """
    Create bw, ip over input
    bwCompare()
    """
    a = read_hiseq(x, hiseq_type)
    if not (a.is_hiseq and hiseq_type.endswith('rx')):
        log.error('hiseq_bw_compare() failed, only support for rx')
        return None
    bw_compare(a.ip_bw, a.input_bw, a.bw, 'subtract',
        threads=a.threads, binsize=100)

