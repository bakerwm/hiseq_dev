#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Arguments for bamCoverage
"""


# import os
# import sys
# import pathlib
import argparse
# import shutil
# import logging
from matplotlib import colors
from hiseq.utils.utils import log


"""
# see: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
# bamCoverage args
bam outFileName outFileFormat MNase Offset region exactScaling
ignoreForNormalization skipNonCoveredRegions ignoreDuplicates
minMappingQuality  samFlagInclude samFlagExclude minFragmentLength maxFragmentLength

scaleFactor binSize normalizeUsing smoothLength extendReads centerReads
blackListFileName effectiveGenomeSize numberOfProcessors filterRNAstrand 
"""


def get_bam_args(**kwargs):
    """
    Arguments for bamCoverage of deeptools
    """
    args = {
        'bam': None,
        'binSize': 50,
        'blackListFileName': None,
        'centerReads': True,
        'effectiveGenomeSize': None,
        'extendReads': None,
        'filterRNAstrand': None,
        'genome': None,
        'normalizeUsing': 'None',
        'numberOfProcessors': 12,
        'out_dir': None,
        'overwrite': False,
        'prefix': 'metaplot',
        'scaleFactor': 1.0,
        'smoothLength': 150,
        'strand_specific': False, # for bigWig files, matrix, ...
    }
    args.update(kwargs)
    return args


def add_common_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return parser    
    parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
        help='the bin_size, default [50]')
    parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
        default=None, help='blacklist file')
    parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4, 
        help='number of processors, default: [4]')
    parser.add_argument('-j', dest='parallel_jobs', type=int, default=1, 
        help='number of jobs run in parallel, default: [4]')
    parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
        help='Overwrite output file')
    # parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
    #     help='labels for samples in plot, default: [None] auto')
    # parser.add_argument('-st', '--startLabel', default='TSS',
    #     help='start label, default: [TSS]')
    # parser.add_argument('-ed', '--endLabel', default='TES',
    #     help='end label, default: [TES]')
    return parser


def add_io_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return parser
    parser.add_argument('-b', dest='bam_list', nargs='+', required=True,
        help='bam files')
    parser.add_argument('-o', dest='out_dir', required=False,
        help='directory to save bigWig file')
    parser.add_argument('-op', '--out-prefix', dest='prefix', default='metaplot',
        help='prefix for output files, default: [metaplot]')
    parser.add_argument('-ss','--strand-specific', dest='strand_specific',
        action='store_true', help='Strand-specific, dUTP library')
    parser = add_common_parser(parser)
    return parser


def add_bam_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return parser
    parser.add_argument('-g', '--genome', default=None,
        help='The reference genome of bam files, default [None]')
    parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
        default=None,
        help='effective genome size, if not specified, parse from bam header')
    parser.add_argument('-s', '--scaleFactor', dest='scaleFactor', type=float,
        default=1.0,
        help='scale factor for the bam, default: [1.0]')
    parser.add_argument('-n', '--normalizeUsing', default='None',
        choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
        help='Use one of the method to normalize reads, default: [None]')
    parser.add_argument('-ex', '--extendReads', type=int, default=None,
        help='extend PE reads to fragment size')
    parser.add_argument('-ce', '--centerReads', action='store_true',
        help='reads are centered with respect to the fragment length')
    parser.add_argument('-sm', '--smoothLength', type=int, default=None,
        help='smooth length')
    parser.add_argument('-fs', '--filterRNAstrand', default=None,
        help='filt RNA strand, forward or reverse')
    return parser

