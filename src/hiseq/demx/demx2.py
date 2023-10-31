#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mission-1. rename fastq files, sub_name to name
Mission-2. demultiplex barcode

## alternative tools
1. defq, ultra fasta multi-threaded fq demultiplexing
https://github.com/OpenGene/defq
2. deML, maxlikelihood demultiplexing
https://github.com/grenaud/deML
"""


import os
import re
import argparse
import tempfile
from difflib import get_close_matches
from pathlib import Path
from xopen import xopen
from collections import Counter
from contextlib import ExitStack # write multiple files
from multiprocessing import Pool
from hiseq.demx.demx_r1 import DemxR1
from hiseq.demx.sample_sheet import HiSeqIndex, SampleSheet
# from sample_sheet import HiSeqIndex, SampleSheet
from hiseq.utils.utils import log, update_obj, Config, get_date, str_distance
from hiseq.utils.seq import Fastx, list_fx, list_fx2, readfq
from hiseq.utils.file import (
    check_dir, file_abspath, file_prefix, file_exists, symlink_file
)


class Demx2(object):
    """
    Rename fastq files, from sub_name to name
    version-1, single i7:      
        # index.csv
        name,i7,i5,bc,reads
        ChrRNA_mHap_DMSO_1h_rep1,TAAGGCGA,NULL,NULL,15
        # raw data
        Next_Ad2.1_1.fq.gz 
        # action
        TAAGGCGA 
            -> Next_Ad2.1 
            -> Next_Ad2.1_1.fq.gz 
            -> ChrRNA_mHap_DMSO_1h_rep1_1.fq.gz (validate i7 index)
    version-2, dual index:
        # index.csv
        sub_name,name,i7,i5,bc,i7_seq,i5_seq,bc_seq,reads
        YY427s001,ChrRNA_mHap_DMSO_1h_rep1,TAAGGCGA,NULL,NULL,NT701,NULL,NULL,15
        # raw data
        YY427s001_1.fq.gz
        # action
        YY427s001_1.fq.gz 
            -> ChrRNA_mHap_DMSO_1h_rep1_1.fq.gz (validate i7+i5 index)

    Arguments
    ---------
    sample_sheet : str
        The sample sheet file, Excel or csv file
        required columns:
        ['sub_name', 'sample_name', 'i7', 'i5', 'bc', 'i7_seq', 'i5_seq',
         'bc_seq', 'reads']
    data_dir, str
        The directory of raw fastq files
    out_dir, str
        The path to the dir, final output
    Example:
    >>> args = {
        'sample_sheet': 'index.csv',
        'data_dir': './from_illumina',
        'out_dir': 'aaa',
    }
    >>> Demx2(**args).run()
    -----------
    This function is designed for renaming fastq files only.
    1. convert `sub_name` to `name`; 
    `sub_name` is a short simplified version for Sequencing Supplier (dry lab)
    `name` is the actual name for the experiment design (wet lab)
    2. validate the i7,i5 index for each fastq file
    3. Statistic the total number of reads, generate a report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'sample_sheet': None,
            'data_dir': None,
            'out_dir': None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'demx_r2'
        try:
            self.data_dir = str(Path(self.data_dir).absolute())
            if not isinstance(self.out_dir, str):
                self.out_dir = Path.cwd()
            self.out_dir = str(Path(self.out_dir).absolute())
        except:
            raise ValueError(f'data_dir not exists: {self.data_dir}')
        self.init_files()
        self.load_fq()
        self.load_index()
        # message
        self.n_fq = len(self.fq_pe_list)
        self.n_idx = len(self.sheet_db)
        print(f'Detected [{self.n_fq}] fastq (pairs), and [{self.n_idx}] index from sample sheet')


    def init_files(self):
        self.sample_sheet = str(Path(self.sample_sheet).absolute())
        self.out_dir = Path(self.out_dir)
        self.config_yaml = str(self.out_dir / 'config.yaml')
        self.index_csv = str(self.out_dir / 'sample_sheet.csv')
        self.fq_metadata_json = str(self.out_dir / 'fq_metadata.json')
        self.fn_json = str(self.out_dir / 'read_count.json')
        self.report_txt = str(self.out_dir / 'report.txt')
        self.out_dir = str(self.out_dir)
        # self.i7_fn_json = self.out_dir / 'i7_index' / 'read_count.json'
        # self.i7_table = self.out_dir / 'index_table.csv'
        # self.bc_table = self.out_dir / 'barcode_table.csv'
        # if not self.config_dir.exists():
        #     self.config_dir.mkdir(parents=True, exist_ok=True)


    def load_fq(self):
        """
        Load fastq files from data_dir
        Expect even number fastq files (paired-end)
        """
        self.fq_info = {} # info in all, see rename_fq(), stat_fq()
        self.fq_list = list_fx(self.data_dir, recursive=True) # sorted !!!
        if len(self.fq_list) > 1: # paired-end mode only
            # split into paired-end mode [[fq1, fq2], [fq1, fq2], ...]
            fq1 = [i for i in self.fq_list if self.fq_paired(i) == 1]
            fq2 = [i for i in self.fq_list if self.fq_paired(i) == 2]
            fq_err = 0
            if len(fq1) == len(fq2):
                for q1, q2 in zip(fq1, fq2):
                    fq_err += str_distance(q1, q2) > 1
            else:
                fq_err = 1 #
            self.fq_is_ok = fq_err == 0
            # pack pe reads
            if self.fq_is_ok:
                self.fq_pe_list = [[i, j] for i,j in zip(fq1, fq2)]
            # guess sheet version:
            # version-1: Next_Ad2.1, TruSeq_Index1, ...
            # version-2: YY427s001, YY427s002, ...
            prefix = [Path(i).stem.lower() for i in self.fq_list] # filename
            pb = [i.startswith('next') or i.startswith('truseq') for i in prefix]
            self.fq_version = 1 if all(pb) else 2 # version 1 or 2
        else:
            # raise ValueError(f'no fastq files found: {self.data_dir}')
            log.error(f'no fastq files found: {self.data_dir}')


    def load_index(self):
        ss = SampleSheet(
            excel_file=self.sample_sheet,
            out_csv=self.index_csv,
            format=self.fq_version, # see self.load_fq()
        )
        ss.to_csv() # convert sheet to csv
        # check if index unique
        dup1 = ss.db.duplicated(subset=['i7', 'i5', 'bc'])
        if any(dup1):
            log.error(ss.db[dup1])
            raise ValueError(f'Duplicate i7/i5/barcode found: {self.sample_sheet}')
        # check if i7 required
        if ss.n_i7 == 0:
            raise ValueError(f'no i7 index found: {self.sample_sheet}')
        # check if i5 exists
        if ss.n_i5 > 0:
            # dual index
            pass
        if ss.n_bc > 0:
            # using demx for each i7/i5 file
            df_bc = ss.db[ss.db['bc'] != 'NULL'] # split 
            raise ValueError(f'use another program for barcode demultiplexing')
        # save as dict
        self.sheet_db = ss.db # pd.DataFrame


    def rename_fq(self):
        """
        rename fastq files (sub_name) by name
        """
        #  check out_dir
        if not Path(self.out_dir).exists():
            Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        if not self.fq_is_ok:
            return None # skipped
        # files in sample_sheet
        ss = self.sheet_db.to_dict('list') # sample_sheet
        sn = ss.get('sub_name', []) # sub_names
        fq_skipped = [] # fq files not in sheet
        for fq1,fq2 in self.fq_pe_list:
            suffix1 = self.fq_suffix(fq1) # _1.fq.gz
            suffix2 = self.fq_suffix(fq2) # _2.fq.gz
            sub_name = self.fq_match(fq1, sn) # sub_name
            if sub_name is None:
                fq_skipped.append(Path(fq1).name)
                continue
            ix = sn.index(sub_name) # index of fq in sheet
            name = ss.get('name')[ix] # new name
            fq1_new = str(Path(self.out_dir) / (name + suffix1))
            fq2_new = str(Path(self.out_dir) / (name + suffix2))
            # print('!A-1', fq1, fq1_new)
            # create symlinks
            # pathlib.relative_to() not working !!!
            if not Path(fq1_new).exists():
                # print('!B-1', fq1_new)
                symlink_file(fq1, fq1_new)
                # Path(fq1_new).symlink_to(fq1) # symlink files
            if not Path(fq2_new).exists():
                # print('!B-2', fq2_new)
                symlink_file(fq2, fq2_new)
                # Path(fq2_new).symlink_to(fq2) # symlink files
            # save info
            self.fq_info.update({
                sub_name : {
                    'sub_name': sub_name,
                    'name': name,
                    'fq1': fq1_new,
                    'fq2': fq2_new,
                    'i7': ss.get('i7')[ix],
                    'i5': ss.get('i5')[ix],
                    'i7_seq': ss.get('i7_seq')[ix],
                    'i5_seq': ss.get('i5_seq')[ix],
                    'reads': ss.get('reads')[ix],
                }
            })
        # warning-1
        if len(fq_skipped) > 0:
            log.warning(f'[{len(fq_skipped)}] fastq files missing from sheet')
            log.warning(','.join(fq_skipped))
        # warning-2
        if len(self.fq_info) != len(sn):
            fq_missing = [i for i in sn if i not in self.fq_info]
            log.warning(f'[{len(fq_missing)}] samples missing fastq files')
            log.warning(','.join(fq_missing))
        # check index and reads
        self.stat_fq()


    def stat_fq(self):
        """
        Wrap all info for each fastq file:
        sub_name, name, index_valid, count
        """
        if not self.fq_is_ok:
            return None # skipped
        # update metadata (from self.fq_info)
        if Path(self.fq_metadata_json).exists():
            self.fq_metadata = Config().load(self.fq_metadata_json)
        else:
            self.fq_metadata = {}
        # show progress
        i = 0
        for k, v in self.fq_info.items(): # see self.rename_fq()
            i += 1
            fq1 = v.get('fq1', None) # fastq file, new_file
            print(f'[{i}/{self.n_fq}] - {fq1}', end='\r') # progress
            # check reads
            fq_count = v.get('count', None)
            if fq_count is None:
                fq_count = Fastx(fq1).number_of_seq()
            # check index
            i7 = v.get('i7_seq', 'NULL')
            i5 = v.get('i5_seq', 'NULL')
            index_valid = v.get('index_valid', None)
            if index_valid is None:
                index_valid = self.validate_index(fq1, i7, i5)
            v.update({
                'index_valid': index_valid,
                'count': fq_count,
            })
            # save as metadata
            self.fq_metadata.update({k:v})
        # save to file
        Config().dump(self.fq_metadata, str(self.fq_metadata_json))


    def fq_match(self, x, p):
        """
        match fastq_name from sub_name
        Get the close match of x from p (list), return the top 1 hit
        x is the query, a string 
        p is the possibilities, a list of sequences
        n is the max number of hits to return, default: 1
        
        Default actions
        Next_Ad2.1: Next_Ad2.1, Next_Ad2-1, Next_Ad2_1
        TruSeq_Index1: TruSeq_Index1

        Output: sub_name
        """
        if isinstance(x, str):
            # pattern: _R1.fastq.gz
            pt = re.compile('(\\.|_)(R?[12]).(f(ast)?q+)(.gz)?$', flags=re.IGNORECASE)
            # prefix
            x_prefix = pt.sub('', str(Path(x).name)) # remove suffix
            hit = get_close_matches(x_prefix, p, n=1, cutoff=0.6) #
            if len(hit) == 0:
                return None
            else:
                return hit[0]
        else:
            return None # no hits


    def fq_paired(self, x):
        """
        Check if fq is paired
        1=read1(PE), 2=read2(PE), 0=SE, -1=(not fastq)
        """
        if isinstance(x, str):
            p1 = re.compile('\\.(f(ast)?q+)(\\.gz)?$', flags=re.IGNORECASE)
            p2 = re.compile('(\\.|_)(R?[12])\\.(f(ast)?q+)(\\.gz)?$', flags=re.IGNORECASE)
            if p1.search(x):
                if p2.search(x):
                    # paired-end
                    r12 = p2.search(x).groups()[1] # R1/2
                    r12 = re.sub('[^0-9]', '', r12) # 1/2
                    r12 = int(r12) #
                else:
                    r12 = 0 # single-end
            else:
                r12 = -1 # not fastq
            return r12


    def fq_suffix(self, x, ext='fq'):
        """
        Get the extension of fastq file
        
        output:
        demo_1.fq.gz     -> _1.fq.gz
        demo_R1.FQ.GZ    -> _1.fq.gz
        demo_R1.FASTQ.gz -> _1.fq.gz
        demo.fq.gz       -> .fq.gz
        demo.fastq.gz    -> .fq.gz
        """
        if isinstance(x, str):
            # suffix: .fq.gz
            if Path(x).suffix == '.gz':
                if isinstance(ext, str):
                    suffix = '.' + ext + '.gz'
                else:
                    suffix = Path(Path(p).stem).suffix + '.gz'
            else:
                if isinstance(ext, str):
                    suffix = '.' + ext
                else:
                    suffix = Path(p).suffix
            # read1/2: _1
            r12 = self.fq_paired(x) # 1=read1, 2=read2, 0=se, -1="not fq"
            if r12 > 0:
                r12 = f'_{r12}' # _1, _2
            else:
                r12 = '' # empty
            # pattern: _R1.fastq.gz
            # pt = re.compile('(\\.|_)(R?[12]).(f(ast)?q+)(.gz)?$', flags=re.IGNORECASE)
            # if pt.search(x):
            #     r12 = pt.search(x).groups()[1] # R1
            #     r12 = re.sub('[^0-9]', '', r12)
            #     r12 = '_' + r12 # _1, _2
            # else:
            #     r12 = ''
            return r12 + suffix


    def validate_index(self, fq, i7, i5, mm=1, n_max=100, cutoff=0.95):
        """
        validate the index for fastq

        Fastq file name_line: (line=1) 
        @A01994:99:HFJCLDSX7:3:1101:19045:1031 1:N:0:GGACTCCT+AGATCTCG
        - 1:N:0:GGACTCCT+AGATCTCG (i7+i5)
        - 1:N:0:ATCACGAT (i7)
        -  (empty) does not contain comment in read_id

        expect 8-digit i5,i7 index
        for 6-digit TruSeq index, append 'AT'

        i7 and i5 is the index_seq
        """
        # hi7 = HiSeqIndex(i7)
        # hi5 = HiSeqIndex(i5)
        # if hi7.is_valid(i7) and hi5.is_valid(i5):
        if isinstance(i7, str) and isinstance(i5, str):
            n_flag = 0
            n_error = 0
            try:
                with xopen(fq) as fh:
                    for r1 in readfq(fh): # name, seq, quality, comment
                        n_flag += 1
                        if n_flag > n_max:
                            break # stop
                        fq_index = list(r1)[-1] # 1:N:0:GGACTCCT+AGATCTCG
                        s = re.split('[:+]', fq_index) # 3:i7, 4:i5
                        m1 = str_distance(s[3], i7) # i7
                        if len(s) > 4:
                            m2 = str_distance(s[4], i5) # i5
                        else:
                            m2 = 0 # single index in fastq
                        if i5 == 'NULL': # skip i5
                            m2 = 0
                        # check mismatch
                        if m1 > mm or m2 > mm:
                            n_error += 1 # error
            except:
                # log.error(f'Could not valid index: {fq}')
                log.warning(f'Could not find index from fq: {fq}')
                n_error = n_max # skipped
            out = (n_max - n_error) / n_max >= cutoff
        else:
            out = False
        return out


    def report(self):
        """
        Generate summary report, save read count
        """
        # save read count to json
        rc = {v.get('name', 'NULL'):v.get('count', 0) for k,v in self.fq_metadata.items()}
        Config().dump(rc, self.fn_json)
        # message
        readsm = sum([v.get('reads', 0) for k,v in self.fq_metadata.items()])
        fq_count = sum([v.get('count', 0) for k,v in self.fq_metadata.items()])
        fq_countm = fq_count / 1e6
        try:
            fq_pct = fq_countm / readsm * 100
        except:
            fq_pct = 100 # divided by zero errors
        # summary
        msg = [
            '='*80,
            f'{"Program":>20} : Demx v2',
            f'{"Date":>20} : {get_date()}',
            f'{"Number of fq(s)":>20} : {self.n_fq}',
            f'{"Number of index":>20} : {self.n_idx}',
            f'{"Total reads":>20} : {fq_count:>11,}',
            f'{"Total reads (M)":>20} : {fq_countm:>11.0f} M',
            f'{"Expect reads (M)":>20} : {readsm:>11} M',
            f'{"Output Percent":>20} : {fq_pct:>11.1f} %',
        ]
        # add header
        msg.append(' '.join([
            f'{"#":>3}',
            f'{"filename":<40}',
            f'{"subname":<16}',
            f'{"i7":<12}',
            f'{"i5":<12}',
            f'{"count":>12}',
            f'{"count(M)":>10}',
            f'{"Exp(M)":>8}',
            f'{"Pct":>5}%',
            f'{"index":>6}'
        ]))
        # metadata
        i = 0
        for k, v in self.fq_metadata.items():
            i += 1
            name = v.get('name', 'NULL')
            subname = v.get('sub_name', 'NULL')
            i7 = v.get('i7', 'NULL')
            i5 = v.get('i5', 'NULL')
            readsm = v.get('reads', 1)
            fq_count = v.get('count', 0)
            fq_countm = fq_count / 1e6
            pct = fq_countm / readsm * 1e2 # %
            index_valid = v.get('index_valid', False)
            index_flag = 'ok' if index_valid else 'failed'
            msg.append(' '.join([
                f'{i:>3}',
                f'{name:<40}',
                f'{subname:<16}',
                f'{i7:<12}',
                f'{i5:<12}',
                f'{fq_count:>12,}',
                f'{fq_countm:>7.1f}   ',
                f'{readsm:>5,}   ',
                f'{pct:>5.1f}',
                f'{index_flag:>6}',
            ]))
        msg.append('='*80)
        # connect
        msg_str = '\n'.join(msg)
        # save msg to file
        with open(self.report_txt, 'wt') as w:
            w.write(msg_str+'\n')
        # to stdout
        print(msg_str)


    def run(self):
        self.rename_fq() # rename files
        self.report()


def get_args():
    example = '\n'.join([
        'De-multiplex version 2',
        '1. rename fastq files, sub_name to name',
        '$ python demx2.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t i7',
        '2. demultiplex barcode',
        '$ python demx2.py -1 fq1.fq.gz -2 fq2.fq.gz -s index.csv -o results -t barcode -x 2 -l 3 -r 2',
    ])
    parser = argparse.ArgumentParser(
        prog='demx2',
        description='De-multiplex version 2',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--data-dir', dest='data_dir', required=True,
        help='Directory contain the fastq files')
    parser.add_argument('-s', '--sample-sheet', dest='sample_sheet', required=True,
        help='sample table in xlsx format, eg: YY00.xlsx')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    # parser.add_argument('-x', '--barcode-in-read', dest='barcode_in_read',
    #     choices=[1, 2], default=2, type=int,
    #     help='barcode in read1/2, default: [2]')
    # parser.add_argument('-l', '--barcode-n-left', type=int,
    #     dest='barcode_n_left', default=0,
    #     help='bases locate on the left of barcode')
    # parser.add_argument('-r', '--barcode-n-right', type=int,
    #     dest='barcode_n_right', default=0,
    #     help='bases locate on the right of barcode')
    # parser.add_argument('-m', '--mismatch', type=int, default=0,
    #     help='mismatches allowed to search index, default: [0]')
    # parser.add_argument('-O', '--overwrite', action='store_true',
    #     help='Overwrite exists files, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx2(**args).run()


if __name__ == '__main__':
    main()

#
