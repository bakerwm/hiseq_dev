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
# from hiseq.demx.sample_sheet import HiSeqIndex, SampleSheet
from sample_sheet import HiSeqIndex, SampleSheet
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
        self.init_files()
        self.load_index()
        self.load_fq()
        # out_dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd())
        # # check fastq files
        # self.raw_fq_list = list_fx(self.data_dir, recursive=True)
        # if len(self.raw_fq_list) < 2:
        #     raise ValueError('no fastq files: {}'.format(self.data_dir))
        # # index table
        # self.idx = self.load_index(self.index_table, min_cols=4)
        # self.demx_type, self.index_list = self.split_i7_bc()
        # if self.demx_type > 1:
        #     self.i7_with_bc = file_prefix(self.index_list[1])
        #     self.i7_with_bc = [i.replace('_table', '') for i in self.i7_with_bc]
        # else:
        #     self.i7_with_bc = [] # empty
        # # save config
        # Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd())
        self.sample_sheet = Path(self.sample_sheet).absolute()
        self.out_dir = Path(self.out_dir).absolute()
        self.config_dir = self.out_dir / 'config'
        self.config_yaml = self.config_dir / 'config.yaml'
        self.index_csv = self.out_dir / 'sample_sheet.csv'
        self.fq_metadata_json = self.out_dir / 'fq_metadata.json'
        self.fn_json = self.out_dir / 'read_count.json'
        # self.i7_fn_json = self.out_dir / 'i7_index' / 'read_count.json'
        self.report_txt = self.out_dir / 'report.txt'
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
        if not self.fq_is_ok:
            return None # skipped
        # files in sample_sheet
        ss = self.sheet_db.to_dcit('list') # sample_sheet
        sn = ss.get('sub_name', []) # sub_names
        fq_skipped = [] # fq files not in sheet
        for fq1,fq2 in self.fq_pe_list:
            suffix1 = self.fq_suffix(fq1) # _1.fq.gz
            suffix2 = self.fq_suffix(fq1) # _2.fq.gz
            sub_name = self.fq_match(fq1, sn) # sub_name
            if sub_name is None:
                fq_skipped.append(Path(fq).name)
                continue
            ix = sn.index(sub_name) # index of fq in sheet
            name = ss.get('name')[ix] # new name
            fq1_new = self.out_dir / name + suffix1
            fq2_new = self.out_dir / name + suffix2
            # create symlinks
            fq1_new.symlink_to(fq1) # symlink files
            fq2_new.symlink_to(fq2) # symlink files
            # save info
            self.fq_info.update({
                sub_name : {
                    'sub_name': sub_name,
                    'name': name,
                    'fq1': fq1_new,
                    'fq2': fq2_new,
                    'i7': ss.get('i7')[ix],
                    'i5': ss.get('i5')[ix],
                    'i7_seq': fq_i7,
                    'i5_seq': fq_i5,
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
        for k, v in self.fq_info.items(): # see self.rename_fq()
            fq1 = v.get('fq1', None) # fastq file, new_file
            # check reads
            fq_count = v.get('count', None)
            if fq_count is None:
                fq_count = Fastx(fq).number_of_seq()
            # check index
            index_valid = v.get('index_valid', None)
            if index_valid is None:
                index_valid = validate_index(fq, i7_seq, i5_seq)
            i7 = v.get('i7_seq', 'NULL')
            i5 = v.get('i5_seq', 'NULL')
            v.update({
                'index_valid': index_valid,
                'count': fq_count,
            })
            # save as metadata
            self.fq_metadata.update({k:v})
        # save to file
        Config().dump(self.fq_metadata, self.fq_metadata_json)


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
                    r12 = p.search(x).groups()[1] # R1/2
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
                        # query = s[3] # i7, ignore i5
                        # m1 = str_distance(s[3], hi7.index) # i7
                        m1 = str_distance(s[3], i7) # i7
                        if len(s) > 4:
                            # m2 = str_distance(s[4], hi5.index) # i5
                            m2 = str_distance(s[4], i5) # i5
                        else:
                            m2 = 0 # single index in fastq
                        # if hi5.index == 'NULL':
                        if i5 == 'NULL': # skip i5
                            m2 = 0
                        # check mismatch
                        if m1 > mm or m2 > mm:
                            n_error += 1 # error
            except:
                log.error(f'Could not read file: {fq}')
                n_error = n_max # skipped
            out = (n_max - n_error) / n_max >= cutoff
        else:
            out = False
        return out


    def report(self):
        """
        Generate summary report, save read count
        """
        # summary
        msg = [
            f'{"Program":>20} : Demx v2',
            f'{"Date":>20} : {get_date()}',
            f'{"Expect reads":>20} : {fq_totalm} M',
            f'{"Total reads":>20} : {fq_totalm:>6.1f} M ({fq_totalm:>11,}) {fq_pct:>8.1f}%',
        ]
        # add header
        msg.append(' '.join([
            f'{"#":>3}',
            f'{"filename":>40}',
            f'{"count":>12}',
            f'{"count(M)":>9}',
            f'{"Expect(M)":>9}',
            f'{"Pct":>5}',
            f'{"index":>5}'
        ])

            #     msg = '\n'.join([
    #         '='*80,
    #         '{:>20} : {}'.format('Program', 'Demx (report)'),
    #         '{:>20} : {}'.format('Date', get_date()),
    #         '{:>20} : {:6.1f} M'.format('Expect reads', total_exp),
    #         '{:>20} : {:6.1f} M ({:>11,}) {:>8.1f}%'.format(
    #             'Total reads', total/1e6, total, output_pct
    #         ),
    #         '{:>5} {:<50s} {:>12} {:>8} {:>8} {:>8}'.format(
    #             'order', 'filename', 'count', 'million', 'expect', 'percent'),
    #         '\n'.join(f_stat),
    #         '='*80,
    #     ])
    #     # save to file
    #     with open(self.report_txt, 'wt') as w:
    #         w.write(msg+'\n')
    #     print(msg)



        i = 0 # number of fastq files
        for k,v in self.fq_metadata.items():
            i += 1
            fq_name = v.get('name', 'NULL')
            fq_reads = v.get('reads', 1) # expect at least 1 M
            fq_count = v.get('count', 0)
            n_pct = fq_count / (fq_reads * 1e4)
            index_valid = v.get('index_valid', False)
            ivalid = 'ok' if index_valid else 'failed'
            # table
            tabs = [
                f'{i:>3}',
                f'{fq_name:<40}',
                f'{fq_count:>12,}',
                f'{fq_count/1e6:>9.1f}',
                f'{fq_reads:>9}',
                f'{n_pct:>5.1f}',
                f'{ivalid:>6},',
            ]
            msg.append(' '.join(tabs))



    def run(self):
        #  check directory
        if not self.config_dir.exists():
            self.config_dir.mkdir(parents=True, exist_ok=True)
        # rename files
        self.rename_fq() # rename files
        # prepare report


    # def report(self):
    #     # load expect read count (million)
    #     df1 = self.load_index(self.index_table, parse_read_count=True)
    #     total_exp = sum(df1.values())
    #     # load real read count
    #     df2 = Config().load(self.fn_json) if file_exists(self.fn_json) else None
    #     if df2 is None:
    #         df2 = {}
    #     total = sum(df2.values())
    #     if total < 1:
    #         total = 1000000 # default: 1M
    #     # check output pct
    #     output_pct = total/(total_exp*1e4) if total_exp > 0 else 100.0
    #     # order
    #     i = 0
    #     f_stat = []
    #     # for k in list(self.idx.values()) + ['undemx']:
    #     idx = self.idx.copy() # local
    #     idx.update({'null':'undemx'})
    #     for k,v in idx.items():
    #         i += 1
    #         n = df2.get(v, 0) # count
    #         n_exp = df1.get(k, 0) # expect read count (million)
    #         n_pct =n/(n_exp*1e4) if n_exp > 0 else 100.0
    #         # elements:
    #         s = ' '.join([
    #             '{:>5d}'.format(i),
    #             '{:<50}'.format(v),
    #             '{:>12,}'.format(n),
    #             '{:>8.1f}'.format(n/1e6),
    #             '{:>8}'.format(n_exp),
    #             '{:>8.1f}%'.format(n_pct),
    #             # '{:8.1f}'.format(v/1e6*scale),
    #         ])
    #         f_stat.append(s)
    #     # output
    #     msg = '\n'.join([
    #         '='*80,
    #         '{:>20} : {}'.format('Program', 'Demx (report)'),
    #         '{:>20} : {}'.format('Date', get_date()),
    #         '{:>20} : {:6.1f} M'.format('Expect reads', total_exp),
    #         '{:>20} : {:6.1f} M ({:>11,}) {:>8.1f}%'.format(
    #             'Total reads', total/1e6, total, output_pct
    #         ),
    #         '{:>5} {:<50s} {:>12} {:>8} {:>8} {:>8}'.format(
    #             'order', 'filename', 'count', 'million', 'expect', 'percent'),
    #         '\n'.join(f_stat),
    #         '='*80,
    #     ])
    #     # save to file
    #     with open(self.report_txt, 'wt') as w:
    #         w.write(msg+'\n')
    #     print(msg)

    # # deprecated
    # def temp_file(self):
    #     # Create a temporary directory
    #     temp_dir = tempfile.TemporaryDirectory()
    #     # Create a temporary file inside the temporary directory
    #     temp_file = Path(temp_dir.name) / "tempfile.txt"
    #     try:
    #         # Write some data to the temporary file
    #         with temp_file.open("w") as f:
    #             f.write("Hello, temporary file!")
    #             # Do some operations with the temporary file
    #         # ...
    #     finally:
    #         # Clean up: Close the temporary directory, which will remove all its contents
    #         temp_dir.cleanup()


    # def load_index(self, x, min_cols=2, parse_read_count=False):
    #     """
    #     index table: name, i7_seq, i5_seq, barcode_seq, ...
    #     ignore: i5 # !!!
    #     read_count (Million), at last column
    #     """
    #     d = {} # i7+bc
    #     try:
    #         with open(x) as r:
    #             for l in r:
    #                 if l.startswith('#') or len(l.strip()) == 0:
    #                     continue
    #                 s = re.split('[,\s\t]', l.strip()) #
    #                 if len(s) < min_cols:
    #                     raise ValueError(
    #                         'at least {} cols required: {}'.format(min_cols, x))
    #                 # match format
    #                 p1 = re.compile('^null$|^[ACGTN]+$', flags=re.IGNORECASE)
    #                 if len(s) >= 4:
    #                     name, i7, i5, bc = s[:4]
    #                     idx = '{}:{}'.format(i7, bc)
    #                     idx_list = s[2:4]
    #                 elif len(s) >= 2:
    #                     name, idx = s[:2]
    #                     idx_list = [s[1]]
    #                 else:
    #                     raise ValueError(
    #                         'at least {} cols required: {}'.format(2, x))
    #                 # if len(s) <= 2:
    #                 #     name, idx = s[:2]
    #                 #     idx_list = [s[1]]
    #                 # elif len(s) <= 4:
    #                 #     name, i7, i5, bc = s[:4]
    #                 #     idx = '{}:{}'.format(i7, bc)
    #                 #     idx_list = s[2:4]
    #                 px = [p1.match(i) is not None for i in idx_list] #
    #                 if not all(px):
    #                     continue
    #                 # parse read count
    #                 p2 = re.compile('^([0-9\.]+)M?$', flags=re.IGNORECASE)
    #                 g2 = p2.match(s[-1])
    #                 if parse_read_count:
    #                     val = eval(g2.group(1)) if g2 else 0
    #                 else:
    #                     val = name
    #                 d.update({idx:val}) # idx:name
    #     except Exception as exc:
    #         log.error(exc)
    #     # check i7
    #     if len(d) == 0:
    #         raise ValueError('no indexes: {}'.format(self.index_table))
    #     return d


    # def split_i7_bc(self):
    #     """
    #     split index table: i7->bc->name
    #     1. into i7
    #     2. into bc
    #     3. into i7 + bc
    #     """
    #     # sub_dir: index table, named by i7
    #     p = re.compile('^null$|^[ACGTN]+$', flags=re.IGNORECASE)
    #     x = [i.split(':') for i in list(self.idx.keys())] # (i7, bc)
    #     i7 = Counter([i[0] for i in x if p.match(i[0])]) # freq
    #     bc = Counter([i[1] for i in x if p.match(i[1])]) # freq
    #     if len(i7) == len(self.idx) or len(bc) == len(self.idx):
    #         tp = 1 # single
    #         idx_list = []
    #         for k,v in self.idx.items():
    #             k1,k2 = k.split(':') # i7,bc
    #             idx_list.append(','.join([v, k1]))
    #         out = self.i7_table if len(i7) == len(self.idx) else self.bc_table
    #         with open(out, 'wt') as w:
    #             w.write('\n'.join(idx_list)+'\n')
    #     else:
    #         # 3. i7 + bc
    #         tp = 2 # multi
    #         ## i7 to N-bc
    #         i7_dir = os.path.join(self.out_dir, 'i7_index')
    #         check_dir(i7_dir)
    #         i7_table = os.path.join(i7_dir, 'i7_table.csv')
    #         out = [i7_table]
    #         i7_list = {} # all unique i7
    #         with open(i7_table, 'wt') as w:
    #             for k,v in self.idx.items():
    #                 k1,k2 = k.split(':') # i7,bc
    #                 if k1 in i7_list:
    #                     continue
    #                 if i7.get(k1, 1) > 1:
    #                     i7v = k1
    #                     i7_list.update({k1:v})
    #                 else:
    #                     i7v = v
    #                 w.write(','.join([i7v, k1])+'\n')
    #         ## N-bc
    #         bc_dir = os.path.join(self.out_dir, 'bc_index')
    #         check_dir(bc_dir) # all
    #         bc_files = [
    #             os.path.join(bc_dir, i+'_table.csv') for i in i7_list
    #         ]
    #         out.append(bc_files)
    #         with ExitStack() as stack:
    #             fws = [stack.enter_context(xopen(f, 'wt')) for f in bc_files]
    #             for k,v in self.idx.items():
    #                 k1,k2 = k.split(':') # i7,bc
    #                 if i7.get(k1, 1) > 1:
    #                     bc_f = os.path.join(bc_dir, k1+'_table.csv')
    #                     fw = fws[bc_files.index(bc_f)]
    #                     fw.write(','.join([v, k2])+'\n')
    #     return [tp, out]


    # def to_i7_index(self, x):
    #     """
    #     index_name: 
    #     TruSeq_Index13
    #     Next_Ad2.1
    #     index_seq:
    #     Index().()
    #     """
    #     # extract i7 name
    #     p = re.compile('True?Seq[\._-]Index\d{1,2}|Next_Ad2[\._-]\d{1,2}|D7\d+', re.IGNORECASE)
    #     g = p.search(x)
    #     if g:
    #         i7 = g.group() # 
    #         i7 = re.sub('Ad2.', 'Ad2.', i7, flags=re.IGNORECASE) # fix Ad2[.-_] to Ad2.
    #         i7 = re.sub('TruSeq.', 'TruSeq_', i7, flags=re.IGNORECASE) # fix TruSeq- to TruSeq_
    #         # Convert i7 seq
    #         out = [i7, HiSeqIndex(i7).index]
    #     else:
    #         out = None
    #     return out


    # def extract_fq_suffix(self, x):
    #     """
    #     Extract the suffix of fastq file:
    #     - *_R1.fq.gz -> _1.fq.gz (PE)
    #     - *_2.fq.gz  -> _2.fq.gz (PE)
    #     - *_r1.fastq -> _1.fq (PE)
    #     """
    #     # xname = os.path.basename(x)
    #     # is_r1 = re.search('_(R)?1.f(ast)?q+.gz', xname, re.IGNORECASE)
    #     p = re.compile('(\.|_)(R?[12]).(f(ast)?q+)(.gz)?', re.IGNORECASE)
    #     g = p.search(x) # 'Next_Ad2-1_R1.fq.gz'
    #     if g:
    #         s = g.groups() # ('_', 'R1', 'fq', None, '.gz')
    #         r12 = re.sub('r', '', s[1], flags=re.IGNORECASE)
    #         out = '_{}.{}'.format(r12, s[2])
    #         if s[4]:
    #             out += s[4]
    #     else:
    #         out = None
    #     return out


    # def rename_by_i7(self, x):
    #     """
    #     fastq file name format:
    #     - Next_Ad2.1_1.fq.gz
    #     - Next_Ad2-1_1.fq.gz
    #     - Next_Ad2_1_1.fq.gz
    #     - YY130-G1-Next_Ad2_1.1.fq.gz
    #     - ...
    #     """
    #     # 1. single i7 mode
    #     i7_df = self.load_index(x, min_cols=2) # i7 table
    #     i7_dir = os.path.join(self.out_dir, 'i7_index')
    #     i7_fn = os.path.join(i7_dir, 'read_count.json')        
    #     i7_fn_df = Config().load(i7_fn) if file_exists(i7_fn) else None
    #     if i7_fn_df is None:
    #         i7_fn_df = {}
    #     check_dir(i7_dir)
    #     # rename files
    #     for fq in self.raw_fq_list:
    #         i7 = self.to_i7_index(fq)
    #         if i7 is None:
    #             # log.error('unknown fq: {}'.format(fq))
    #             continue
    #         i7_id, i7_seq = i7
    #         name = i7_df.get(i7_seq, None)
    #         suffix = self.extract_fq_suffix(fq)
    #         if name is None:
    #             # log.error('unknown i7: {} {}'.format(i7_id, i7_seq))
    #             continue
    #         if suffix is None:
    #             # log.error('unknown fq suffix: {}'.format(fq))
    #             continue
    #         new_fq = os.path.join(i7_dir, name+suffix)
    #         if not file_exists(new_fq):
    #             symlink_file(fq, new_fq)
    #         # count reads
    #         fq_count = i7_fn_df.get(name, -1)
    #         if fq_count < 0:
    #             # fq_count = Fastx(new_fq).number_of_seq()
    #             if i7_seq not in self.i7_with_bc:
    #                 fq_count = Fastx(new_fq).number_of_seq()
    #             else:
    #                 fq_count = 1
    #             i7_fn_df.update({name:fq_count})
    #         if i7_seq not in self.i7_with_bc:
    #             if not suffix.startswith('_2'): # skip read2
    #                 log.info('check file: {} {}'.format(name, fq_count))
    #     # update read_count.json
    #     Config().dump(i7_fn_df, i7_fn)


    # def demx_bc(self, x):
    #     """
    #     Parameters:
    #     ----------
    #     x : str
    #         barcode index table
    #     """
    #     args = self.__dict__.copy()
    #     i7_dir = os.path.join(self.out_dir, 'i7_index')
    #     i7_fq_list = list_fx(i7_dir, recursive=False)
    #     bc = file_prefix(x).replace('_table', '') # barcode seq
    #     bc_fq = [i for i in i7_fq_list if os.path.basename(i).startswith(bc)]
    #     args.update({
    #         'fq1': bc_fq[0] if len(bc_fq) > 0 else None,
    #         'fq2': bc_fq[1] if len(bc_fq) > 1 else None,
    #         'out_dir': os.path.join(self.out_dir, 'bc_index', bc),
    #         'index_type': 'barcode',
    #         'index_table': x,
    #     })
    #     DemxR1(**args).run() #


    # def demx_i7_bc(self):
    #     args = self.__dict__.copy()
    #     tp, tb = self.split_i7_bc()
    #     if tp == 1:
    #         # 1 run i7
    #         self.rename_by_i7(tb)
    #     else:
    #         # 2. multi i7/bc mode
    #         tb_i7, tb_bc = tb #
    #         ## 2.1 run i7
    #         self.rename_by_i7(tb_i7)
    #         ## 2.2 run bc
    #         n_bc = 8 if len(tb_bc) > 8 else len(tb_bc)
    #         if n_bc > 0:
    #             with Pool(processes=n_bc) as pool:
    #                 pool.map(self.demx_bc, tb_bc)


    # def wrap_i7_bc(self):
    #     tp, tb = self.split_i7_bc()
    #     # i7
    #     i7_df = Config().load(self.i7_fn_json) if file_exists(self.i7_fn_json) else None
    #     if i7_df is None:
    #         i7_df = {}
    #     i7_undemx = i7_df.get('undemx', 0)
    #     # barcode
    #     bc_df = {}
    #     i7_drop = []
    #     if tp > 1:
    #         for bc_f in tb[1]:
    #             bc = file_prefix(bc_f)
    #             bc = bc.replace('_table', '') # barcode seq
    #             bc_fc = os.path.join(self.out_dir, 'bc_index', bc, 'read_count.json')
    #             df = Config().load(bc_fc) #
    #             i7_undemx += df.get('undemx', 0)
    #             bc_df.update(df) # barcode files
    #             [i7_drop.append(k) for k,v in i7_df.items() if k.startswith(bc)]
    #     # merge
    #     i7_df.update(bc_df)
    #     [i7_df.pop(i, None) for i in i7_drop]
    #     i7_df.update({'undemx': i7_undemx}) # update undemx
    #     Config().dump(i7_df, self.fn_json)
    #     # i7+bc: rename files
    #     i7_dir = os.path.join(self.out_dir, 'i7_index')
    #     bc_dir = os.path.join(self.out_dir, 'bc_index')
    #     for i in list(self.idx.values()):
    #         # i7 index
    #         i7_list = list_fx2(i7_dir, i+'*', recursive=False)
    #         # bc index
    #         if file_exists(bc_dir):
    #             bc_list = list_fx2(bc_dir, i+'*', recursive=True)
    #         else:
    #             bc_list = []
    #         for a in i7_list + bc_list:
    #             a_new = os.path.join(self.out_dir, os.path.basename(a))
    #             if not file_exists(a_new):
    #                 symlink_file(a, a_new)


    # def show_msg(self):
    #     msg = '\n'.join([
    #         '='*80,
    #         '{:>20} : {}'.format('Program', 'Demx_index'),
    #         '{:>20} : {}'.format('Date', get_date()),
    #         '{:>20} : {}'.format('fq_dir', self.data_dir),
    #         '{:>20} : {}'.format('index_table', self.index_table),
    #         '{:>20} : {}'.format('out_dir', self.out_dir),
    #         '{:>20} : {}'.format('report', self.report_txt),
    #         '{:>20} : {}'.format('mismatch', self.mismatch),
    #         '{:>20} : {}'.format('overwrite', 'yes' if self.overwrite else 'no'),
    #         '='*80,
    #     ])
    #     print(msg)


    # def report(self):
    #     # load expect read count (million)
    #     df1 = self.load_index(self.index_table, parse_read_count=True)
    #     total_exp = sum(df1.values())
    #     # load real read count
    #     df2 = Config().load(self.fn_json) if file_exists(self.fn_json) else None
    #     if df2 is None:
    #         df2 = {}
    #     total = sum(df2.values())
    #     if total < 1:
    #         total = 1000000 # default: 1M
    #     # check output pct
    #     output_pct = total/(total_exp*1e4) if total_exp > 0 else 100.0
    #     # order
    #     i = 0
    #     f_stat = []
    #     # for k in list(self.idx.values()) + ['undemx']:
    #     idx = self.idx.copy() # local
    #     idx.update({'null':'undemx'})
    #     for k,v in idx.items():
    #         i += 1
    #         n = df2.get(v, 0) # count
    #         n_exp = df1.get(k, 0) # expect read count (million)
    #         n_pct =n/(n_exp*1e4) if n_exp > 0 else 100.0
    #         # elements:
    #         s = ' '.join([
    #             '{:>5d}'.format(i),
    #             '{:<50}'.format(v),
    #             '{:>12,}'.format(n),
    #             '{:>8.1f}'.format(n/1e6),
    #             '{:>8}'.format(n_exp),
    #             '{:>8.1f}%'.format(n_pct),
    #             # '{:8.1f}'.format(v/1e6*scale),
    #         ])
    #         f_stat.append(s)
    #     # output
    #     msg = '\n'.join([
    #         '='*80,
    #         '{:>20} : {}'.format('Program', 'Demx (report)'),
    #         '{:>20} : {}'.format('Date', get_date()),
    #         '{:>20} : {:6.1f} M'.format('Expect reads', total_exp),
    #         '{:>20} : {:6.1f} M ({:>11,}) {:>8.1f}%'.format(
    #             'Total reads', total/1e6, total, output_pct
    #         ),
    #         '{:>5} {:<50s} {:>12} {:>8} {:>8} {:>8}'.format(
    #             'order', 'filename', 'count', 'million', 'expect', 'percent'),
    #         '\n'.join(f_stat),
    #         '='*80,
    #     ])
    #     # save to file
    #     with open(self.report_txt, 'wt') as w:
    #         w.write(msg+'\n')
    #     print(msg)


    # def run(self):
    #     self.show_msg()
    #     self.demx_i7_bc()
    #     self.wrap_i7_bc()
    #     self.report()


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
    parser.add_argument('-s', '--index-table', dest='index_table',
        required=True,
        help='index table in csv format, [filename,barcode]')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    parser.add_argument('-x', '--barcode-in-read', dest='barcode_in_read',
        choices=[1, 2], default=2, type=int,
        help='barcode in read1/2, default: [2]')
    parser.add_argument('-l', '--barcode-n-left', type=int,
        dest='barcode_n_left', default=0,
        help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int,
        dest='barcode_n_right', default=0,
        help='bases locate on the right of barcode')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def main():
    args = vars(get_args().parse_args())
    Demx2(**args).run()


if __name__ == '__main__':
    main()

#
