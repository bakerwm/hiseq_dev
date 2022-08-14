#!/usr/bin/env python
#-*- encoding:utf-8 -*-

"""
Convert bam to bigwig: bamCoverage
# strand-specific:
fwd: scale(fwd): fwd / total
rev: scale(rev): rev / total
# output
outdir/out.bigWig
outdir/out_fwd.bigWig
outdir/out_rev.bigWig
Example:
$ bamCoverage -b in.bam -o out.bw --binSize 50 \
  --effectiveGenomeSize 100 \
  --scaleFactor 1.0 --normalizeUsing None \
  --blackListFileName b.bed \
  --skipNAs --extendReads --centerReads -p 8
"""


import os
import pathlib
import argparse
import shutil
import pysam
from multiprocessing import Pool
from hiseq.utils.utils import log, update_obj, Config, init_cpu
from hiseq.utils.genome import Genome
from hiseq.bam2bw.bam2bw_args import get_bam_args, add_io_parser, add_bam_parser
from hiseq.utils.file import (
    file_abspath, file_prefix, file_exists, symlink_file, fix_out_dir, 
    is_valid_file, is_valid_bam, is_valid_bigwig
)


class Bam2bw(object):
    def __init__(self, **kwargs):
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args()


    def update_args(self):
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list] # to list
        if not is_valid_file(self.bam_list, is_valid_bam):
            msg = '\n'.join([
                '{} : {}'.format(is_valid_file(i, is_valid_bam), i) for i in self.bam_list
            ])
            print(msg)
            raise ValueError('bam file illegal')


    def run_single_bam(self, i):
        func = Bam2bw_ss if self.strand_specific else Bam2bw_ns
        args = self.__dict__.copy()
        args.update({'bam': self.bam_list[i], 'prefix': None, })
        return func(**args).run()

        
    def run(self):
        # run in parallel
        if len(self.bam_list) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                out = pool.map(self.run_single_bam, range(len(self.bam_list)))
        else:
            out = [self.run_single_bam(i) for i in range(len(self.bam_list))]
        return out # bw list
        

class Bam2bw_ns(object):
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        # self.update_labels()
        self.init_files()
        self.cmd = self.get_cmd()
        Config().dump(self.__dict__, self.config)


    def basic_args(self):
        """
        default arguments [38]
        to-do: outFileSortedRegions, outFileNameMatrix
        """
        # bam, outFileName, outFileFormat
        arg_list = [
            'bam', 'outFileName',
            'scaleFactor', 'normalizeUsing', 'binSize', 'effectiveGenomeSize',
            'blackListFileName', 'numberOfProcessors',
            'extendReads', 'smoothLength',
            'filterRNAstrand', 'region', 'MNase', 'Offset', 
            'exactScaling', 'ignoreForNormalization',
            'minMappingQuality', 'samFlagInclude', 
            'samFlagExclude', 'minFragmentLength', 'maxFragmentLength'
        ]
        # ['skipNAs', 'centerReads', 'ignoreDuplicates'] # no arguments
        return arg_list


    def update_args(self):
        arg_list = self.basic_args()
        d = {i:getattr(self, i, None) for i in arg_list}
        self = update_obj(self, d, force=True) # update
        if not isinstance(self.prefix, str):
            self.prefix = file_prefix(self.bam)
        if not is_valid_file(self.bam, is_valid_bam):
            raise ValueError('bam file illegal: {}'.format(self.bam))
        # effsize
        es = getattr(self, 'effectiveGenomeSize', None)
        if es is None:
            es = self.get_effsize()
            setattr(self, 'effectiveGenomeSize', es)
        # blacklist
        bl = getattr(self, 'blackListFileName', None)
        if isinstance(self.genome, str):
            if bl is None:
                self.blackListFileName = Genome(self.genome).blacklist()


    def get_effsize(self):
        """
        1. get effective genomesize from default values
        2. return the chromsize from bam file
        """
        effsize = {
            'dm3': 162367812,
            'dm6': 142573017,
            'mm9': 2620345972,
            'mm10': 2652783500,
            'hg19': 2451960000,
            'hg38': 2913022398,
            'GRCh38': 2913022398
        }
        # genome
        genome = getattr(self, 'genome', None)
        s = effsize.get(genome, None)
        # from bam file
        if not isinstance(s, int):
            sam = pysam.AlignmentFile(self.bam)
            s = sum(sam.header.lengths)
        if not isinstance(s, int):
            raise ValueError('could not determine effsize: {}'.format(self.bam))
        return s


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.prefix)
        if not file_exists(self.project_dir):
            os.makedirs(self.project_dir)
        self.bam = file_abspath(self.bam)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'bw': prefix+'.bigWig',
            'stdout': prefix+'.bamCoverage.stdout',
            'stderr': prefix+'.bamCoverage.stderr',
            'bw_cmd': prefix+'.bamCoverage.sh',
        }
        self = update_obj(self, args, force=True)
        self.outFileName = self.bw #


    def get_cmd(self):
        """
        construct arguments to command line
        """
        alist = self.basic_args()
        args = {i:getattr(self, i, None) for i in alist}
        dlist = ['--{} {}'.format(k, v) for k,v in args.items() if v is not None]
        bb = ['skipNAs', 'centerReads', 'ignoreDuplicates'] # no arguments
        bba = ['--'+i for i in bb if getattr(self, i, None)]
        dlist += bba # add arguments
        dline = ' '.join(dlist) # to cmd line
        # main args
        cmd = ' '.join([
            '{}'.format(shutil.which('bamCoverage')),
            dline,            
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        return cmd


    def run(self):
        with open(self.bw_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # bam index
        bai = self.bam + '.bai'
        if not os.path.exists(bai):
            pysam.index(self.bam)
        # run
        if os.path.exists(self.bw) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('bamCoverage() skipped, file exists: {}'.format(self.bw))
        else:
            log.info('run bamCoverage: {}'.format(self.bw))
            os.system(self.cmd)
        # check output
        if not is_valid_file(self.bw, is_valid_bigwig):
            log.error('bamCoverage() failed, file not found: {}'.format(self.bw))
        return self.bw


class Bam2bw_ss(object):
    """
    Strand-specific:
    add scale: fwd (fwd/fwd+rev); rev (rev/fwd+rev)
    """
    def __init__(self, **kwargs):
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        self.init_files()
        Config().dump(self.__dict__, self.config)


    def update_args(self):
        prefix = getattr(self, 'prefix', None)
        if not isinstance(prefix, str):
            self.prefix = file_prefix(self.bam)
        if not is_valid_file(self.bam, is_valid_bam):
            raise ValueError('bam file illegal: {}'.format(self.bam))


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.prefix)
        self.project_dir = fix_out_dir(self.project_dir)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'bw': prefix+'.bigWig',
            'bw_fwd': prefix+'_fwd.bigWig',
            'bw_rev': prefix+'_rev.bigWig',
        }
        self = update_obj(self, args, force=True)


    def count_bam(self, bam, strand=None):
        """
        for strand-specific RNA-seq (NSR)
        read2 is the sense direction
        -f 16 : forward
        -F 16 : reverse
        """
        c1 = 0 
        c2 = 0
        if strand == '+' or strand is None:
            # (forward) -–filterRNAstrand=forward keeps minus-strand reads, -f 16
            c1 = pysam.view('-c', '-f', '16', '-F', '4', '-@', '8', bam) # fwd
            c1 = int(c1.strip())
        elif strand == '-' or strand is None:
            # (reverse) -–filterRNAstrand=reverse keeps plus-strand reads, -F 16
            c2 = pysam.view('-c', '-F', '16', '-F', '4', '-@', '8', bam) # rev
            c2 = int(c2.strip())
        else:
            log.error('unknown strand: {}'.format(strand))
        return c1 + c2


    def get_bam_count(self, bam, strand=None):
        """
        get BAM count from file, or samtools view -c
        save count to yaml: count.yaml
        """
        bname = file_prefix(bam)
        cf = os.path.join(self.project_dir, bname+'.count.yaml')
        if not os.path.exists(cf):
            cfwd = self.count_bam(bam, strand='+')
            crev = self.count_bam(bam, strand='-')
            d = {
                'name': bname,
                'forward': cfwd,
                'reverse': crev,
                'total': cfwd + crev,
            }
            Config().dump(d, cf)
        # get from yaml
        d = Config().load(cf)
        dn = d.get('name', None)
        if strand == '+':
            out = d.get('forward', 0)
        elif strand == '-':
            out = d.get('reverse', 0)
        else:
            out = d.get('total', 0)
        return out


    def get_bam_scale(self, bam):
        """
        the norm scale for forward/reverse strand
        """
        fwd = self.get_bam_count(self.bam, strand='+')
        rev = self.get_bam_count(self.bam, strand='-')
        sf = round(fwd/(fwd+rev), 6)
        sr = round(rev/(fwd+rev), 6)
        return (sf, sr)


    def run(self):
        sf, sr = self.get_bam_scale(self.bam)
        # bam index
        bai = self.bam + '.bai'
        if not os.path.exists(bai):
            pysam.index(self.bam)
        # forward
        args1 = self.__dict__.copy()
        args1.update({
            'out_dir': os.path.join(self.out_dir, self.prefix),
            'prefix': self.prefix+'_fwd',
            'scaleFactor': sf,
            'filterRNAstrand': 'forward',
            # 'normalizeUsing': 'CPM',
        })
        b1 = Bam2bw_ns(**args1)
        b1.run()
        # reverse
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': os.path.join(self.out_dir, self.prefix),
            'prefix': self.prefix+'_rev',
            'scaleFactor': sf,
            'filterRNAstrand': 'reverse',
            # 'normalizeUsing': 'CPM',
        })
        b2 = Bam2bw_ns(**args2)
        b2.run()
        # link files to upper-level folder
        symlink_file(b1.bw, self.bw_fwd)
        symlink_file(b2.bw, self.bw_rev)
        # check output
        if not is_valid_file([self.bw_fwd, self.bw_rev], is_valid_bigwig):
            log.error('bamCoverage() failed, {} {}'.format(self.bw_fwd, self.bw_rev))

        return (self.bw_fwd, self.bw_rev)


def get_args():
    example = ' '.join([
        '$ bamCoverage -b in.bam -o out.bw --binSize 50',
        '--effectiveGenomeSize 100',
        '--scaleFactor 1.0 --normalizeUsing None',
        '--blackListFileName b.bed',
        '--skipNAs --extendReads 150 --centerReads -p 8',
    ])
    parser = argparse.ArgumentParser(
        prog='bam2bw.py', description='bam2bw', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = add_io_parser(parser)
    parser = add_bam_parser(parser)
    return parser


def main():
    args = vars(get_args().parse_args())
    Bam2bw(**args).run()


if __name__ == '__main__':
    main()

