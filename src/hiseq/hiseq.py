#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is the main script for hiseq
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08 at hotmail.com'
__date__ = '2021-05-20'
__version__ = '1.0.1'


import os
import sys
import argparse
from hiseq.utils.utils import log

## functions
from hiseq.download.download import download
from hiseq.download.download import get_args as add_download_args

from hiseq.demx.sample_sheet import SampleSheet
from hiseq.demx.sample_sheet import get_args as add_sheet_args

from hiseq.demx.demx import Demx
from hiseq.demx.demx import get_args as add_demx_args

from hiseq.demx.demx2 import Demx2
from hiseq.demx.demx2 import get_args as add_demx2_args

from hiseq.trim.trim import Trim
from hiseq.trim.trim import get_args as add_trim_args


# from hiseq.qc.fastqc import Fastqc
# from hiseq.qc.fastqc import get_args as add_fastqc_args

# from hiseq.align.align import Align
# from hiseq.align.align import get_args as add_align_args

# from hiseq.qc.hiseq_lib import HiseqLib, HiseqLibSmRNA
# from hiseq.qc.hiseq_lib import get_args as add_p7_args

# # sub-modules-1
# from hiseq.bam2bw.bam2bw import Bam2bwRn
# from hiseq.bam2bw.bam2bw import get_args as add_bam2bw_args

# from hiseq.utils.bam import Bam2cor, Bw2cor
# from hiseq.bam2cor.bam2cor import get_args as add_bam2cor_args

# from hiseq.fragsize.fragsize import BamFragSize
# from hiseq.fragsize.fragsize import get_args as add_fragsize_args

# from hiseq.utils.hiseq_list import hiseq_list, list_dirs
# from hiseq.utils.hiseq_list import get_args as add_hiseq_list_args

# # from hiseq.qc.parse_i7 import HiSeqP7
# # from hiseq.qc.parse_i7 import get_args as add_p7_args
# from hiseq.qc.bacteria import Kraken2
# from hiseq.qc.bacteria import get_args as add_bacteria_args

# from hiseq.get_trackhub.get_trackhub import TrackHub
# from hiseq.get_trackhub.get_trackhub import get_args as add_trackhub_args



# # to-be-deprecated: replaced by specific get_args() in each command
# from hiseq.utils.argsParser import add_quant_args, add_peak_args, add_motif_args, \
#     add_rnaseq_args2, add_deseq_pair_args, add_go_args,\
#     add_peak2idr_args, add_bed2overlap_args, \
#     add_sample_args

# from hiseq.utils.fastx import Fastx
# from hiseq.sample.sample import FxSample
# from hiseq.utils.utils import log, get_date, print_dict
# from hiseq.utils import download as dl


class Hiseq(object):
    """
    The 1st-level of command, choose which sub-command to use
    qc, trim, align, quant, peak, motif, report, ... (to be continued)
    """
    def __init__(self):
        usage = """ hiseq <command> [<args>]

    subcommands:

        sheet        Preparing sample_sheet.csv for demx/demx2
        demx         Demultiplexing reads (P7, barcode) from single Lane
        demx2        Demultiplexing multi barcode files
        qc           quality control, fastqc
        p7           Check the P7 of HiSeq library
        bacteria     Check bacteria content

        atac         ATACseq pipeline
        cnr          CUN&RUN pipeline
        chipseq      ChIPseq pipeline
        rnaseq       RNAseq pipeline
        rnaseq2      RNAseq pipeline, simplify version
        rnaseq_salmon RNAseq pipeline using Salmon+DESeq2

        trim         trim adapters, low-quality bases, ...
        trim_smRNA   trim adapters, UMI, for small RNA libraries
        align        Align fastq/a files to reference genome
        quant        Count genes/features
        peak         Call peaks using MACS2
        motif        Check motifs from a BED/fasta file
        report       Create a report to the above commands
        go           Run GO analysis on geneset
        deseq_pair   Run RNAseq compare

        run_trackhub Generate the trackhub urls
        fragsize     Fragment size of PE alignments
        bam2cor      Correlation between bam files
        bam2bw       Convert bam to bigWig
        peak2idr     Calculate IDR for multiple Peaks
        bed2overlap  Calculate the overlap between bed intervals
        sample       Sample fastq file
        download     Download files
    """
        parser = argparse.ArgumentParser(
            prog = 'hiseq',
            description = 'A collection of tools for HiSeq data',
            epilog = '',
            usage = usage,
        )
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            log.error('unknown command: {}'.format(args.command))
            parser.print_help()
            sys.exit(1)
        getattr(self, args.command)()


    def init_args(self, p):
        """
        Parameters
        ----------
        p:
            argparse.ArgumentPaser
        """
        if isinstance(p, argparse.ArgumentParser):
            out = vars(p.parse_args(sys.argv[2:]))
        else:
            out = None
        return out


## sub-modules
    def download(self):
        args = self.init_args(add_download_args())
        download(**args)


    def sheet(self):
        """
        Prepare the sample sheet for demx
        1. sample_name,i7,i5,barcode (Demx, whole-lane)
        2. i7_name,i7,reads (Demx2, part)
        """
        args = self.init_args(add_sheet_args())
        SampleSheet(**args).run()


    def demx(self):
        args = self.init_args(add_demx_args())
        Demx(**args).run()


    def demx2(self):
        args = self.init_args(add_demx2_args())
        Demx2(**args).run()


    def trim(self):
        args = self.init_args(add_trim_args())
        Trim(**args).run()




# ## pipelines

#     def atac(self):
#         """
#         ATACseq pipeline
#         """
#         args = self.init_args(add_atac_args())
#         Atac(**args).run()


#     def cnr(self):
#         """
#         CUN&RUN pipeline
#         """
#         args = self.init_args(add_cnr_args())
#         Cnr(**args).run()


#     def hiseq_merge(self):
#         """
#         Merge multiple hiseq dirs
#         """
#         args = self.init_args(add_hiseq_merge_args())
#         CnrMerge(**args).run()


#     def chipseq(self):
#         args = self.init_args(add_chipseq_args())
#         Chipseq(**args).run()



#     def rnaseq(self):
#         """
#         RNA-seq pipeline
#         """
#         args = self.init_args(add_rnaseq_args())
#         Rnaseq(**args).run()


#     def rnaseq_salmon(self):
#         """
#         RNA-seq pipeline, using salmon
#         """
#         args = self.init_args(add_rnaseq_salmon_args())
#         RnaseqSalmonPipe(**args).run()


#     def trim_smRNA(self):
#         """
#         small RNA pipeline
#         """
#         args = self.init_args(add_trim_smRNA_args())
#         TrimSmRNA(**args).run()


#     def trim(self):
#         """
#         Trim adapters
#         """
#         args = self.init_args(add_trim_args())
#         Trim(**args).run()


#     def align(self):
#         """
#         Align reads
#         """
#         args = self.init_args(add_align_args())
#         Align(**args).run()


#     def bam2bw(self):
#         """
#         Convert bam to bw files
#         """
#         args = self.init_args(add_bam2bw_args())
#         Bam2bwRn(**args).run()


#     def bam2cor(self):
#         """
#         Calculate bam correlation
#         using deeptools
#         """
#         args = self.init_args(add_bam2cor_args())
# #         Bam2cor(**args).run()
#         if all([i.endswith('.bam') for i in args['bam_list']]):
#             print('Bam2cor')
#             Bam2cor(**args).run()
#         elif all([i.endswith('.bigWig') for i in args['bam_list']]):
#             print('Bw2cor')
#             args['bw_list'] = args['bam_list'] # update
#             Bw2cor(**args).run()
#         else:
#             log.error('no bam/bigWig files found')


#     def bed2overlap(self):
#         """
#         Calculate IDR for peak files
#         using: idr
#         """
#         args = self.init_args(add_bed2overlap_args())
#         BedOverlap(**args).run()


#     def fragsize(self):
#         """
#         Calculate the fragment size of PE alignment
#         """
#         args = self.init_args(add_fragsize_args())
#         BamFragSize(**args).run()




#     def qc(self):
#         """
#         Fastq quality control
#         """
#         args = self.init_args(add_fastqc_args())
#         Fastqc(**args).run()


#     def hiseq_list(self):
#         args = self.init_args(add_hiseq_list_args())
#         dirs = args.pop('dirs', None)
#         dirs = list_dirs(dirs) # update, root/subdir
#         out = hiseq_list(dirs, **args)
#         print('\n'.join(out))


#     def p7(self):
#         """
#         Check library structure: P7, barcode
#         """
#         args = self.init_args(add_p7_args())
# #         HiseqLib(**args).run()
#         if args['smRNA']:
#             HiseqLibSmRNA(**args).run()
#         else:
#             HiseqLib(**args).run()


#     def bacteria(self):
#         """
#         Check bacteria content
#         """
#         args = self.init_args(add_bacteria_args())
#         Kraken2(**args).run()


#     def run_trackhub(self):
#         """
#         Make trackhub
#         """
#         args = self.init_args(add_trackhub_args())
#         TrackHub(**args).run()


# ################################################################################
# ## to-be-updated
# ## 2021-05-20
#     def quant(self):
#         """
#         quantify hiseq reads
#         """
#         args = self.init_args(add_quant_args())
#         print('hiseq quant')


#     def peak(self):
#         """
#         Call peaks
#         """
#         args = self.init_args(add_peak_args())
#         CallPeak(**args).run()


#     def motif(self):
#         """
#         quantify hiseq reads
#         ...
#         """
#         args = self.init_args(add_motif_args())
#         print('hiseq motif')


#     def go(self):
#         """
#         Run GO analysis on geneset
#         """
#         args = self.init_args(add_go_args())
#         if args['all'] is None and args['input'] is None:
#             parser.parse_args(['-h'])
#             sys.exit('arguments failed, either --all or --input required')
#         Go(**args).run()


#     def deseq_pair(self):
#         """
#         Run RNAseq cmp
#         """
#         args = self.init_args(add_deseq_pair_args())
#         DeseqPair(**args).run()


# #     def rnaseq2(self):
# #         """
# #         RNA-seq pipeline, simplify version
# #         """
# #         args = self.init_args(add_rnaseq_args2)
# #         RNAseqPipe(**args).run()


# #     def chipseq(self):
# #         """
# #         ChIPseq pipeline
# #         """
# #         args = self.init_args(add_chipseq_args())
# #         ChIPseq(**args).run()


#     def peak2idr(self):
#         """
#         Calculate IDR for peak files
#         using: idr
#         """
#         args = self.init_args(add_peak2idr_args())
#         PeakIDR(**args).run()


#     def sample(self):
#         """
#         Sub-sample fastq files
#         head
#         """
#         args = self.init_args(add_sample_args())
#         FxSample(**args).run()


def main():
    Hiseq()


if __name__ == '__main__':
    main()

#