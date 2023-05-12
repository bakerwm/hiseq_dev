#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions for BAM files
Bam
# Bam2cor
# Bam2fingerprint
BAM:
  - sort
  - index
  - count
  - rm_dup
  - proper_pair
  - subset
  - frag_size *
  - read_size *
  - to_bed
  - to_bg *
  - to_bw *
  - to_fq *
  - to_xx *
"""


import os
import tempfile
import numpy as np
import pysam
import pybedtools
import shutil
import pathlib
import re
from shutil import which
from hiseq.utils.utils import log, update_obj, run_shell_cmd, Config
from hiseq.utils.file import check_dir, check_file, file_exists, file_abspath, file_prefix, list_file, read_file


class Bam(object):
    """
    Operation on BAM files
    - sort
    - index
    - merge
    - count
    - to_bed
    - rm_dup
    - ...
    Using Pysam,...
    code from cgat:
    """
    def __init__(self, bam, threads=4, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.bam = bam
        self.threads = threads


    def index(self):
        bai = self.bam + '.bai'
        if not os.path.exists(bai):
            pysam.index(self.bam)
        # return os.path.exists(bai)


    def sort(self, outfile=None, by_name=False, overwrite=False):
        """
        Sort bam file by position (default)
        save to *.sorted.bam (or specify the name)
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.sorted.bam'
        if os.path.exists(outfile) and overwrite is False:
            log.info('file exists: {}'.format(outfile))
        else:
            if by_name:
                tmp = pysam.sort('-@', str(self.threads), '-n', '-o', outfile, self.bam)
            else:
                tmp = pysam.sort('-@', str(self.threads), '-o', outfile, self.bam)
        return outfile


    def merge(self):
        """
        Merge multiple BAM files using samtools
        """
        # pysam.merge('')
        pass


    def count(self, reads=True):
        x = pysam.view('-c', self.bam)
        n = int(x.strip())
        if not reads:
            if self.is_paired():
                n = int(n / 2)
        return n


    def to_bed(self, outfile=None):
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.bed'
        if not os.path.exists(outfile):
            pybedtools.BedTool(self.bam).bam_to_bed().saveas(outfile)
        return outfile


    def rm_dup(self, outfile=None, overwrite=False, tools='picard'):
        """
        Remove duplicates using picard/sambamba
        sambamba markdup -r --overflow-list-size 800000 raw.bam rm_dup.bam
        picard MarkDuplicates -REMOVE_DUPLICATES True -I in.bam -O outfile.bam -M metrix.txt
        ## -REMOVE_SEQUENCING_DUPLICATES True
        """
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.rm_dup.bam'
        if tools == 'sambamba':
            sambamba = which('sambamba')
            log_stderr = outfile + '.sambamba.log'
            cmd = ' '.join([
                '{} markdup -r'.format(which('sambamba')),
                '-t {}'.format(2), #!! force to 2
                '--overflow-list-size 1000000',
                '--tmpdir={}'.format(tempfile.TemporaryDirectory().name),
                '{} {} 2> {}'.format(self.bam, outfile, log_stderr),
            ])
        elif tools == 'picard':
            picard = which('picard')
            log_stderr = outfile + '.picard.log'
            metrics_file = outfile + '.metrics.txt'
            cmd = ' '.join([
                '{} MarkDuplicates'.format(which('picard')),
                'REMOVE_DUPLICATES=True',
                'I={} O={} M={}'.format(self.bam, outfile, metrics_file),
                '2>{}'.format(log_stderr),
                '&& samtools index {}'.format(outfile)
            ])
            # 'REMOVE_SEQUENCING_DUPLICATES=True',
        else:
            log.error(' '.join([
                'rm_dup(), unknown tools {}, '.format(tools),
                'options: ["sambamba", "picard"]'
            ]))
            return None
        # save cmd
        cmd_txt = outfile.replace('.bam', '.cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # if os.path.exists(outfile) and overwrite is False:
        if check_file(outfile, check_empty=True) and not overwrite:
            log.info('file exists: {}'.format(outfile))
        else:
            run_shell_cmd(cmd)
        return outfile


    def proper_pair(self, outfile=None, overwrite=False):
        if outfile is None:
            outfile = os.path.splitext(self.bam)[0] + '.proper_pair.bam'
        if os.path.exists(outfile) and overwrite is False:
            log.info('file exists: {}'.format(outfile))
        else:
            pysam.view('-f', '2', '-h', '-b', '-@', str(self.threads),
                '-o', outfile, self.bam, catch_stdout=False)
        return outfile


    def subset(self, size=20000, subdir=None):
        if subdir is None:
            subdir = self._tmp(delete=False)
        check_dir(subdir)
        # src, dest
        dest = os.path.join(subdir, os.path.basename(self.bam))
        if file_exists(dest):
            log.info('Bam.subset() skipped, file exists {}'.format(dest))
        else:
            self.index()
            srcfile = pysam.AlignmentFile(self.bam, 'rb')
            destfile = pysam.AlignmentFile(dest, 'wb', template=srcfile)
            # counter
            i = 0
            for read in srcfile.fetch():
                i +=1
                if i > size:
                    break
                destfile.write(read)
        return dest


    def _tmp(self, delete=True):
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=delete)
        return tmp.name


    ##########################################
    ## code from cgat: BEGIN
    ##########################################
    def is_paired(self, topn=1000):
        """
        Check if bam contains paired end reads
        go through the topn alignments in file,
        return: True, any of the alignments are paired
        """
        samfile = pysam.AlignmentFile(self.bam)
        n = 0
        for read in samfile:
            if read.is_paired:
                break
            n += 1
            if n == topn:
                break
        samfile.close()
        return n != topn


    def get_num_reads(self):
        """
        Count number of reads in bam file.
        This methods works through pysam.idxstats.
        Arguments
        ---------
        bamfile : string
            Filename of :term:`bam` formatted file. The file needs
            to be indexed.
        Returns
        -------
        nreads : int
            Number of reads
        """
        lines = pysam.idxstats(self.bam).splitlines()
        try:
            nreads = sum(
                map(int, [x.split("\t")[2]
                          for x in lines if not x.startswith("#")]))
        except IndexError as msg:
            raise IndexError(
                "can't get number of reads from bamfile, msg=%s, data=%s" %
                (msg, lines))
        return nreads


    def estimateInsertSizeDistribution(self, topn=10000, n=10,
        method="picard", similarity_threshold=1.0, max_chunks=1000):
        """
        Estimate insert size from a subset of alignments in a bam file.
        Several methods are implemented.
        picard
            The method works analogous to picard by restricting the estimates
            to a core distribution. The core distribution is defined as all
            values that lie within n-times the median absolute deviation of
            the full data set.
        convergence
            The method works similar to ``picard``, but continues reading
            `alignments` until the mean and standard deviation stabilize.
            The values returned are the median mean and median standard
            deviation encountered.
        The method `convergence` is suited to RNA-seq data, as insert sizes
        fluctuate significantly depending on the current region
        being looked at.
        Only mapped and proper pairs are considered in the computation.
        Returns
        -------
        mean : float
           Mean of insert sizes.
        stddev : float
           Standard deviation of insert sizes.
        npairs : int
           Number of read pairs used for the estimation
        method : string
           Estimation method
        similarity_threshold : float
           Similarity threshold to apply.
        max_chunks : int
           Maximum number of chunks of size `alignments` to be used
           in the convergence method.
        """
        assert self.is_paired(self.bam), \
            'can only estimate insert size from' \
            'paired bam files'
        samfile = pysam.AlignmentFile(self.bam)
        def get_core_distribution(inserts, n):
            # compute median absolute deviation
            raw_median = np.median(inserts)
            raw_median_dev = np.median(np.absolute(inserts - raw_median))
            # set thresholds
            threshold_min = max(0, raw_median - n * raw_median_dev)
            threshold_max = raw_median + n * raw_median_dev
            # define core distribution
            return inserts[np.logical_and(inserts >= threshold_min,
                                          inserts <= threshold_max)]
        if method == "picard":
            # only get first read in pair to avoid double counting
            inserts = np.array(
                [read.template_length for read in samfile.head(n=topn)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)
            return np.mean(core), np.std(core), len(inserts)
        elif method == "convergence":
            means, stds, counts = [], [], []
            last_mean = 0
            iteration = 0
            while iteration < max_chunks:
                inserts = np.array(
                    [read.template_length for read in samfile.head(
                        n=topn,
                        multiple_iterators=False)
                     if read.is_proper_pair
                     and not read.is_unmapped
                     and not read.mate_is_unmapped
                     and not read.is_read1
                     and not read.is_duplicate
                     and read.template_length > 0])
                core = get_core_distribution(inserts, n)
                means.append(np.mean(core))
                stds.append(np.std(core))
                counts.append(len(inserts))
                mean_core = get_core_distribution(np.array(means), 2)
                mm = np.mean(mean_core)
                if abs(mm - last_mean) < similarity_threshold:
                    break
                last_mean = mm
            return np.median(means), np.median(stds), sum(counts)
        else:
            raise ValueError("unknown method '%s'" % method)


    def estimateTagSize(self, topn=10, multiple="error"):
        """
        Estimate tag/read size from first alignments in file.
        Parameters
        ---------
        bamfile : string
           Filename of :term:`bam` formatted file
        alignments : int
           Number of alignments to inspect
        multiple : string
           How to deal if there are multiple tag sizes present.
           ``error`` will raise a warning, ``mean`` will return the
           mean of the read lengths found. ``uniq`` will return a
           unique list of read sizes found. ``all`` will return all
           read sizes encountered.
        Returns
        -------
        size : int
           The read size (actual, mean or list of read sizes)
        Raises
        ------
        ValueError
           If there are multiple tag sizes present and `multiple` is set to
           `error`.
        """
        samfile = pysam.AlignmentFile(self.bam)
        sizes = [read.rlen for read in samfile.head(topn)]
        mi, ma = min(sizes), max(sizes)
        if mi == 0 and ma == 0:
            sizes = [read.inferred_length for read in samfile.head(topn)]
            # remove 0 sizes (unaligned reads?)
            sizes = [x for x in sizes if x > 0]
            mi, ma = min(sizes), max(sizes)
        if mi != ma:
            if multiple == "error":
                raise ValueError('multiple tag sizes in %s: %s' % (self.bam, sizes))
            elif multiple == "mean":
                mi = int(sum(sizes) / len(sizes))
            elif multiple == "uniq":
                mi = list(sorted(set(sizes)))
            elif multiple == "all":
                return sizes
        return mi


    def getNumberOfAlignments(self):
        """Return number of alignments in bamfile.
        """
        if not os.path.exists(self.bam + '.bai'):
            pysam.index(self.bam)
        samfile = pysam.AlignmentFile(self.bam)
        return samfile.mapped


def is_sam_flag(x, return_codes=False):
    """
    For sam flags
    The flags for sam format are in binary code:
    see: https://samtools.github.io/hts-specs/SAMv1.pdf
    1    0x1   template having multiple segments in sequencing
    2    0x2   each segment properly aligned according to the aligner
    4    0x4   segment unmapped
    8    0x8   next segment in the template unmapped
    16   0x10  SEQ being reverse complemented
    32   0x20  SEQ of the next segment in the template being reverse complemented
    64   0x40  the first segment in the template
    128  0x80  the last segment in the template
    256  0x100 secondary alignment
    512  0x200 not passing filters, such as platform/vendor quality controls
    1024 0x400 PCR or optical duplicate
    2048 0x800 supplementary alignment
    see: Bitwise operator
    """
    f = {
        '1': 'template having multiple segments in sequencing',
        '2': 'each segment properly aligned according to the aligner',
        '4': 'segment unmapped',
        '8': 'next segment in the template unmapped',
        '16': 'SEQ being reverse complemented',
        '32': 'SEQ of the next segment in the template being reverse complemented',
        '64': 'the first segment in the template',
        '128': 'the last segment in the template',
        '256': 'secondary alignment',
        '512': 'not passing filters, such as platform/vendor quality controls',
        '1024': 'PCR or optical duplicate',
        '2048': 'supplementary alignment',
    }
    # for details, valid code
    h = []
    if isinstance(x, int):
        if x in range(1, 2049):
            for k,v in f.items():
                k = int(k)
                if k == k & x:
                    x -= k
                    h.append((k, v))
        else:
            log.error('x not valid, expect [1, 2048], got: {}'.format(x))
    else:
        log.error('x not valid, expect int, got {}'.format(type(x).__name__))
    if return_codes and h:
        out = h
    else:
        out = len(h) > 0
    return out


class Bam2cor(object):
    """
    Compute correlation between replicates
    input: bam
    output: count_matrix
            cor_matrix
            cor_plot
            ...
    window = 500bp
    eg:
    multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
        --outRawCounts *counts.tab -b bam
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'cor_method': 'pearson',
            'make_plot': True,
            'bam_dir': None,
            'bam_list': None,
            'out_dir': None,
            'prefix': None,
            'binsize': 500,
            'threads': 1,
            'overwrite': False,
            'multibamsummary': shutil.which('multiBamSummary'),
            'plotcorrelation': shutil.which('plotCorrelation'),
            'plotpca': shutil.which('plotPCA'),
            'flag': True # whether run/not
        }
        self = update_obj(self, args_init, force=False)
        # output
        if self.out_dir is None:
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        if not isinstance(self.prefix, str):
            self.prefix = 'multibam'
        # check commands
        self.flag = all([
            isinstance(self.multibamsummary, str),
            isinstance(self.plotcorrelation, str),
            isinstance(self.plotpca, str),
        ])
        if not self.flag:
            raise ValueError('check bam and commands')
        self.init_bam()
        self.init_files()
        # save config
        check_dir(self.out_dir, create_dirs=True)
        Config().dump(self.__dict__, self.config_yaml)


    def init_bam(self):
        b = []
        if isinstance(self.bam_dir, str):
            b = list_file(self.bam_dir, '*.bam')
        elif isinstance(self.bam_list, str):
            if self.bam_list.endswith('.bam'):
                b.append(self.bam_list)
            elif self.bam_list.endswith('.txt'):
                b = read_file(self.bam_list, comment='#')
            else:
                log.error('unknown bam: {}'.format(self.bam_list))
        elif isinstance(self.bam_list, list):
            b = self.bam_list
        else:
            log.error('bam_dir, bam required')
        # file exists
        b = [i for i in b if file_exists(i) and i.endswith('.bam')]
        [Bam(i).index() for i in b] #
        if len(b) == 0:
            log.error('no bam files detected')
        self.bam_list = b
        self.bam_line = ' '.join(self.bam_list)
        [Bam(b).index() for b in self.bam_list]


    def init_files(self):
        default_files = {
            'config_yaml': os.path.join(self.out_dir, 'config.yaml'),
            'bam_npz': os.path.join(self.out_dir, self.prefix + '.npz'),
            'log': os.path.join(self.out_dir, self.prefix + '.deeptools.log'),
            'plot_cor_heatmap_png': os.path.join(self.out_dir, self.prefix + '.cor_heatmap.png'),
            'log_heatmap': os.path.join(self.out_dir, self.prefix + '.cor_heatmap.log'),
            'cor_counts': os.path.join(self.out_dir, self.prefix + '.cor_counts.tab'),
            'cor_matrix': os.path.join(self.out_dir, self.prefix + '.cor.matrix'),
            'plot_cor_pca_png': os.path.join(self.out_dir, self.prefix + '.cor_PCA.png'),
            'log_pca': os.path.join(self.out_dir, self.prefix + '.cor_PCA.log')
        }
        self = update_obj(self, default_files, force=True) # key
        check_dir(self.out_dir, create_dirs=True)


    def bam_summary(self):
        cmd = ' '.join([
            '{} bins --binSize {}'.format(self.multibamsummary, self.binsize),
            '-p {}'.format(self.threads),
            '--smartLabels -o {}'.format(self.bam_npz),
            '--outRawCounts {}'.format(self.cor_counts),
            '--bamfiles {}'.format(self.bam_line)
        ])
        cmd_txt = os.path.join(self.out_dir, self.prefix + '.bam_cor.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(self.bam_npz) and not self.overwrite:
            log.warning('file exists: {}'.format(self.bam_npz))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
        # check
        if not os.path.exists(self.bam_npz):
            log.error('Bam2cor() failed, output file not found: {}'.format(
                self.bam_npz))


    def plot_cor_heatmap(self):
        cmd = ' '.join([
            '{} -in {}'.format(self.plotcorrelation, self.bam_npz),
            '--corMethod {}'.format(self.cor_method),
            '--plotTitle "{} Correlation"'.format(self.cor_method),
            '--whatToPlot heatmap --colorMap RdYlBu --plotNumbers',
            '-o {}'.format(self.plot_cor_heatmap_png),
            '--outFileCorMatrix {}'.format(self.cor_matrix)
            ])
        cmd_txt = os.path.join(self.out_dir, self.prefix + '.heatmap.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(self.plot_cor_heatmap_png) and not self.overwrite:
            log.warning('Bam2cor.plot_cor_heatmap() skipped, file exists: {}'.format(
                self.plot_cor_heatmap_png))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log_heatmap, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
        # check
        if not os.path.exists(self.plot_cor_heatmap_png):
            log.error('Bam2cor() failed, output file not found: {}'.format(
                self.plot_cor_heatmap_png))


    def plot_cor_pca(self):
        """
        plotPCA:
        Make PCA plot, for bam_summary
        correlation
        """
        cmd = ' '.join([
            '{} -in {}'.format(self.plotpca, self.bam_npz),
            '-o {}'.format(self.plot_cor_pca_png),
            '-T "PCA for Bam files"'
            ])
        cmd_txt = os.path.join(self.out_dir, self.prefix + '.pca.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(self.plot_cor_pca_png) and not self.overwrite:
            log.warning('Bam2cor.plot_cor_pca() skipped, file exists: {}'.format(self.bam_npz))
        else:
            _, stdout, stderr = run_shell_cmd(cmd)
            with open(self.log_pca, 'wt') as w:
                w.write(stdout + '\n' + stderr + '\n')
        # check
        if not os.path.exists(self.plot_cor_pca_png):
            log.error('Bam2cor() failed, output file not found: {}'.format(self.plot_cor_pca_png))


    def run(self):
        if len(self.bam_list) > 1:
            self.bam_summary()
            # plot
            if self.make_plot is True:
                self.plot_cor_heatmap()
                self.plot_cor_pca()
        else:
            log.error('require >=2 bam files')


class Bam2fingerprint(object):
    """
    QC for bam files
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        Required arguments for ATACseq analysis
        """
        args_init = {
            'bam_dir': None,
            'bam_list': None,
            'threads': 8,
            'out_dir': None,
            'prefix': None,
            'title': 'BAM_fingerprint',
            'labels': None,
            'overlap': False,
            'overwrite': False
            }
        self = update_obj(self, args_init, force=False)
        if self.out_dir is None:
            self.out_dir = str(pathlib.Path.cwd())
        self.init_bam()
        self.init_files()
        # save config
        check_dir(self.out_dir, create_dirs=True)
        Config().dump(self.__dict__, self.config_yaml)


    def init_bam(self):
        b = []
        if isinstance(self.bam_dir, str):
            b = list_file(self.bam_dir, '*.bam')
        elif isinstance(self.bam_list, str):
            if self.bam_list.endswith('.bam'):
                b.append(self.bam_list)
            elif self.bam_list.endswith('.txt'):
                b = read_file(self.bam_list, comment='#')
            else:
                log.error('unknown bam: {}'.format(self.bam_list))
        elif isinstance(self.bam_list, list):
            b = self.bam_list
        else:
            log.error('bam_dir, bam required')
        # file exists
        b = [i for i in b if file_exists(i) and i.endswith('.bam')]
        [Bam(i).index() for i in b] #
        if len(b) == 0:
            log.error('no bam files detected')
        self.bam_list = b
        self.bam_line = ' '.join(self.bam_list)
        [Bam(b).index() for b in self.bam_list]


    def init_files(self):
        if not (isinstance(self.labels, list) and len(self.labels) == len(self.bam_list)):
            self.labels = file_prefix(self.bam_list)
        self.labels_line = ' '.join(self.labels)
        ## output files
        prefix = re.sub('[^A-Za-z0-9-.]', '_', self.title)
        prefix = re.sub('_+', '_', prefix) # format string
        if isinstance(self.prefix, str):
            prefix = self.prefix
        else:
            self.prefix = prefix
        self.fp_png = os.path.join(self.out_dir, prefix + '.png')
        self.fp_tab = os.path.join(self.out_dir, prefix + '.tab')
        self.config_yaml = os.path.join(self.out_dir, 'config.yaml')


    def run_fingerprint(self):
        cmd = ' '.join([
            '{}'.format(shutil.which('plotFingerprint')),
            '--minMappingQuality 30 --skipZeros --numberOfSamples 50000',
            '-p {}'.format(self.threads),
            '--plotFile {}'.format(self.fp_png),
            '--outRawCounts {}'.format(self.fp_tab),
            '-b {}'.format(self.bam_line),
            '--labels {}'.format(self.labels_line)
            ])
        # save cmd
        cmd_txt = os.path.join(self.out_dir, self.prefix + '_cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # run
        if file_exists(self.fp_png) and not self.overwrite:
            log.info('Bam2fingerprint() skipped, file exists: {}'.format(self.fp_png))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('Bam2fingerprint() failed, see {}'.format(self.out_dir))


    def run(self):
        msg = '\n'.join([
            '-'*80,
            '{:>14} : {}'.format('program', 'hiseq.utils.bam.Bam2fingerprint()'),
            '{:>14} : {}'.format('bam', self.bam_list),
            '{:>14} : {}'.format('n_bam', len(self.bam_list)),
            '{:>14} : {}'.format('out_dir', self.out_dir),
            '{:>14} : {}'.format('prefix', self.prefix),
            '-'*80,
        ])
        print(msg)
        self.run_fingerprint()
                
