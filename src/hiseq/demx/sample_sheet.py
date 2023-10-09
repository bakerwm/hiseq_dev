#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The header of sheet:
file: sample_sheet.xlsx 
version: 2023-10-08
1.SampleID
2.Lib_number
3.Lib_sub
4.Lib_user
5.Sample_name
6.RBP
7.Cell_line
8.Species
9.Spike-in
10.P7_index_id
11.P5_index_id
12.Barcode_id
13.Seq_type
14.Lib_type
15.Reads,M
16.Reads,G (FCID)
17.sub_name (Lane)
18.Reference
19.Spikein_ref
20.P7_index
21.P5_index
22.Barcode_seq
23.Reminder
24.Status
25.Date

version: 2022-04-26
1.SampleID
2.Lib_number
3.Lib_sub
4.Lib_user
5.Sample_name
6.RBP
7.Cell_line
8.Species
9.Spike-in
10.P7_index_id
11.Barcode_id
12.Seq_type
13.Lib_type
14.Reads,M
15.FCID
16.Lane
17.Reference
18.Spikein_ref
19.P7_index
20.Barcode_seq
21.Reminder
22.Status
23.Date

Version: 2021-06-01
header:
1.SampleID
2.Lib_number
3.Lib_user
4.Sample_name
5.RBP
6.Cell_line
7.Species
8.Spike-in
9.P7_index_id
10.Barcode_id
11.Seq_type
12.Lib_type
13.Reads,M
14.FCID
15.Lane
16.Reference
17.Spikein_ref
18.P7_index
19.Barcode_seq

# input: sheet.xlsx
# output: demx.csv
sub_name,name,i7,i5,bc,i7_seq,i5_seq,bc_seq,reads
YY278s001,CnT_hela_pcf11_sc514158_rep1,Next_Ad2.1,NULL,NULL,TAAGGCGA,NULL,NULL,15

# sub-functions
1. fix_sample_name()
2. index_to_seq()
3. seq_to_index()

Mission:
1. Sample_name, sanitize, unique
2. P7_index_id, Barcode_id, unique
3. index_table: p7,p5,bc
4. index_table: p7,bc (sub_table)
5. p7 to name
6. reads(M)
to_MGI:
p7_id,p7_seq,reads(Gb)
"""


import os
import sys
import re
# import hiseq
# import pathlib
import warnings
from pathlib import Path
import argparse
import pandas as pd
# from hiseq.utils.helper import update_obj
from hiseq.utils.utils import log, update_obj, get_date
import importlib.resources as res

# mission-1
# 1. convert to MGI table
# 2. split by i7 index


class HiSeqIndex(object):
    """
    Convert between index_name and index_seq
    Unique: name to seq
    multi-hits: seq to name (version, lib_type)
    Example:
    >>> from sample_sheet import HiSeqIndex
    ## input index_seq
    >>> i = HiSeqIndex('ATCACG')
    >>> i.name
    'TruSeq_Index1'
    >>> i.index
    'ATCACG'
    ## input index_name
    >>> i = HiSeqIndex('TruSeq_Index1')
    >>> i.name
    'TruSeq_Index1'
    >>> i.index
    'ATCACG'
    ## input invalid
    """
    def __init__(self, x):
        self.x = x
        self.db = self.load_hiseq_index() # name:{"id": id, "seq": seq, ...}
        # self.df2 = {v:k for k,v in self.df1.items()} # seq:name
        if isinstance(x, str):
            self.name = self.seq_to_name(x)
            self.index = self.name_to_seq(x)
        elif isinstance(x, list):
            self.name = [self.seq_to_name(i) for i in x]
            self.index = [self.name_to_seq(i) for i in x]
        else:
            raise ValueError(f'illegal x={x}, expect str,list got {type(x)}')


    def load_hiseq_index(self):
        # locate the index file
        with res.path('hiseq.data', 'illumina_index.csv') as f:
            ff = str(f)
        d = {}
        with open(ff) as r:
            for line in r:
                if line.startswith('#'):
                    continue
                tabs = line.strip().split(',')
                if len(tabs) < 6:
                    continue
                id, name, seq, idx, lib, version = tabs[:6]
                if id == 'index_id': # header line
                    continue
                # fix version
                version = re.sub('[^0-9.]', '', version) # numeric
                try:
                    version = float(version) # convert to numeric
                except:
                    version = 1.0 # default
                d.update({
                    id:{
                        "id": id,
                        "name": name,
                        "seq": seq,
                        "idx": idx,
                        "lib": lib,
                        "version": version
                    }
                })
        return d


    def is_valid(self, x):
        if isinstance(x, str):
            out = self.is_valid_name(x) or self.is_valid_index(x)
        elif isinstance(x, list):
            out = [self.is_valid(i) for i in x]
        else:
            out = None # unknown
        return out


    def is_valid_name(self, x):
        return x in self.db or str(x).upper() == 'NULL'


    def is_valid_index(self, x):
        s = [v.get('seq', None) for k,v in self.db.items()] # seq list
        return x in s or str(x).upper() == 'NULL'


    def seq_to_name(self, x, version='latest', lib='auto'):
        # lib: TruSeq > NSR
        if self.is_valid_name(x):
            out = x
        elif self.is_valid_index(x):
            dx = [] # candidates
            for k,v in self.db.items():
                v_seq = v.get('seq', None)
                v_version = v.get('version', None)
                v_lib = v.get('lib', None)
                # check if seq
                if x == v_seq:
                    dx.append(v) #
            # check if multiple hits
            if len(dx) == 1:
                out = dx[0].get('id', 'NULL')
            else: # more than 1 hits
                # check version
                ver_list = [i.get('version', None) for i in dx]
                if version not in ver_list:
                    version = max(ver_list)
                # check lib: TruSeq > NSR
                lib_list = [i.get('version', None) for i in dx]
                if lib not in lib_list:
                    lib = 'TruSeq'
                # check version and lib
                di = [i.get('id', 'NULL') for i in dx 
                      if i.get('version', 0) == version 
                      and i.get('lib', None) == lib]
                out = di[0]
        else:
            out = 'NULL'
        return out


    def name_to_seq(self, x):
        if self.is_valid_index(x):
            out = x
        elif self.is_valid_name(x):
            df = self.db.get(x, None)
            out = df.get('seq', None)
        else:
            out = 'NULL'
        return out


class SampleSheet(object):
    """
    Parsing data from Excel sample sheet file.xlsx
    # Example
    >>> s = SampleSheet(excel_file='YY00.xlsx', out_csv='demx.csv')
    >>> s.to_barcode_table()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'excel_file': None,
            'out_csv': None,
            'format': 2,
        }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.excel_file, str):
            raise ValueError(f'-i not str, got {type(self.excel_file)}')
        if not Path(self.excel_file).exists():
            raise ValueError(f'-i not exists: {self.excel_file}')
        if not Path(self.excel_file).suffix in ['.xlsx', '.xls', '.csv']:
            raise ValueError(f'-i not [xlsx, xls, csv]: {self.excel_file}')
        if not isinstance(self.out_csv, str):
            raise ValueError(f'-o not str, got {type(self.out_csv)}')
        self.excel_file = Path(self.excel_file).absolute() # absolute path
        self.out_csv = Path(self.out_csv).absolute() # absolute path
        self.out_dir = str(Path(self.out_csv).parent)
        self.db = self.read_table(self.excel_file) #
        self.stat_table()


    def read_table(self, x):
        """
        Sample sheet format:
        version-1:
        sample_name, p7_index_id, p5_index_id, barcode_id, readsm
        version-2:
        sub_name, sample_name, p7_index_id, p5_index_id, barcode_id, readsm
        output: pd.DataFrame
        """
        try:
            if Path(x).suffix in ['.csv']:
                df = pd.read_csv(x, comment='#')
            elif Path(x).suffix in ['.xlsx', 'xls']:
                # suppress warnings of Data validation by openpyxl
                warnings.simplefilter(action='ignore', category=UserWarning) 
                df = pd.read_excel(x, sheet_name='sample_sheet')
                warnings.resetwarnings() # reset to default
            else:
                pass
            # update table
            df.dropna(axis='index', thresh=4, inplace=True) # require 4 columns
            df.fillna('NULL', inplace=True) # convert nan to 'NULL'
            # df.replace('^Null$', 'NULL', regex=True, inplace=True)
            # require regex, ignore case ?
            df.replace(
                regex=['^Null$', '^NUll$', '^NULl$', '^NUlL$', '^NuLL$'], 
                value='NULL', inplace=True
            )
            c0 = df.columns.to_list() # original columns
            df.columns = [re.sub('[^\\w]', '', i).lower() for i in c0]
        except:
            log.error(f'Could not read sample_sheet: {x}')
            df = pd.DataFrame(
                columns=['sub_name', 'name', 'i7', 'i5', 'bc', 'reads']
            )
        return self.fix_table(df)


    def fix_table(self, df):
        """
        format table, get index name, index seq, sample name
        """         
        # headers
        # version-1:
        header1 = ['sample_name', 'p7_index_id', 'barcode_id', 'readsm']
        header1b = ['name', 'i7', 'bc', 'reads']
        # version-2:
        header2 = [
            'sub_name', 'sample_name', 'p7_index_id', 'p5_index_id',
            'barcode_id', 'readsm'
        ]
        header2b = ['sub_name', 'name', 'i7', 'i5', 'bc', 'reads']
        # final
        header3b = [
            'sub_name', 'name', 'i7', 'i5', 'bc', 'i7_seq', 'i5_seq',
            'bc_seq', 'reads'
        ]
        # check: version-2
        if all([i in df.columns for i in header3b]):
            return df # final format
        if all([i in df.columns for i in header2b]):
            pass
        elif all([i in df.columns for i in header2]):
            df = df[header2] # subset
            df.columns = header2b # rename header
        # check: version-1
        elif all([i in df.columns for i in header1b]):
            df = df[header1b] # subset
            df.loc[:, ['sub_name']] = df['name']
            df.loc[:, ['i5']] = 'NULL'
        elif all([i in df.columns for i in header1]):
            df = df[header1] # subset
            df.columns = header1b # rename header
            df.loc[:, ['sub_name']] = df['name']
            df.loc[:, ['i5']] = 'NULL'
        else:
            raise ValueError(f'unkown sample sheet format: {x}')
        # format:
        if all(df['name'] == df['sub_name']) or self.format == 1:
            df.loc[:, ['sub_name']] = df['i7'] # i7_name for sub_name
        # fix sample name format
        df.loc[:, ['name']] = self.sanitize(df['name'].to_list())
        # index_name to index_seq
        df = df.assign(
            i7_seq = HiSeqIndex(df['i7'].to_list()).index,
            i5_seq = HiSeqIndex(df['i5'].to_list()).index,
            bc_seq = HiSeqIndex(df['bc'].to_list()).index
        )
        return df[header3b] # arrange columns


    def stat_table(self):
        """
        Statistics of the excel table
        """        
        # sample name
        self.name_is_unique = sum(self.db.duplicated(subset=['name'])) == 0
        self.sub_name_is_unique = sum(self.db.duplicated(subset=['sub_name'])) == 0
        # total sample
        self.n_sample = self.db.shape[0]
        # Number of i7
        db_7 = self.db[self.db['i7'] != 'NULL']
        self.n_i7 = db_7.groupby('i7').ngroups
        # Number of i5
        db_5 = self.db[self.db['i5'] != 'NULL']
        self.n_i5 = db_5.groupby('i5').ngroups
        # Number of barcode
        db_bc = self.db[self.db['bc'] != 'NULL']
        self.n_bc = db_bc.groupby('bc').ngroups
        # Number of reads, million
        self.n_reads = self.db['reads'].sum()

    
    # deprecated: see read_table()
    def read_csv(self, x):
        """
        Sample sheet in csv format:
        version-1:
        sample_name, p7_index_id, p5_index_id, barcode_id, readsm
        version-2:
        sub_name, sample_name, p7_index_id, p5_index_id, barcode_id, readsm
        output: pd.DataFrame
        """
        try:
            df = pd.read_csv(x, comment='#')
            # df.dropna(axis='index', thresh=4, inplace=True) # require 4 columns
            # df.fillna('NULL', inplace=True) # convert nan to 'NULL'
            # #(?i) ignore case
            # df.replace('^(?i)null$', 'NULL', regex=True, inplace=True)
            # c0 = df.columns.to_list() # original columns
            # df.columns = [re.sub('[^\w]', '', i).lower() for i in c0]
        except:
            df = pd.DataFrame(
                columns=[
                    'sub_name', 'sample_name', 'p7_index_id', 'p5_index_id',
                    'barcode_id', 'readsm'
                ]
            )
        return df

    
    # deprecated: see read_table()
    def read_excel(self, x):
        try:
            # suppress warnings of Data validation by openpyxl
            warnings.simplefilter(action='ignore', category=UserWarning) 
            df = pd.read_excel(x, sheet_name='sample_sheet')
            warnings.resetwarnings() # reset to default
            # df.dropna(axis='index', thresh=4, inplace=True) # require 4 columns
            # df.fillna('NULL', inplace=True) # convert nan to 'NULL'
            # #(?i) ignore case
            # df.replace('^(?i)null$', 'NULL', regex=True, inplace=True)
            # c0 = df.columns.to_list() # original columns
            # df.columns = [re.sub('[^\w]', '', i).lower() for i in c0]
        except:
            df = pd.DataFrame(
                columns=[
                    'sub_name', 'sample_name', 'p7_index_id', 'p5_index_id',
                    'barcode_id', 'readsm'
                ]
            )
        return df


    def sanitize(self, x):
        """
        fix sample_name, remove non alphanumeric characters
        1. [A-Za-z0-9_-.]
        2. convert to '_'
        """
        if isinstance(x, str):
            out = re.sub('[^\\w\\-.]', '', x)
            out = re.sub('_+', '_', out)
        elif isinstance(x, list):
            out = [self.sanitize(i) for i in x]
        else:
            log.warning('illegal x')
            out = x
        return out


    def to_csv(self):
        """
        save as csv file, with following columns
        sub_name,name,i7,i5,bc,i7_seq,i5_seq,bc_seq,reads
        """
        try:
            if not Path(self.out_dir).exists():
                Path(self.out_dir).mkdir(parents=True, exist_ok=True)
            self.db.to_csv(self.out_csv, index=False, header=True)
        except:
            log.error(f'Could not write to file: {self.out_csv}')


    def run(self):
        # db = self.read_table(self.excel_file) #
        self.to_csv() # save to csv file
        msg = '\n'.join([
            '='*80,
            f'{"Program":<20}: {"SampleSheet":}',
            f'{"Date":<20}: {get_date():}',
            f'{"Excel file":<20}: {self.excel_file:}',
            f'{"Demx csv file":<20}: {self.out_csv:}',
            f'{"Output format":<20}: version-{self.format}',
            f'{"No. of samples":<20}: {self.n_sample:}',
            f'{"No. of i7":<20}: {self.n_i7:}',
            f'{"No. of i5":<20}: {self.n_i5:}',
            f'{"No. of barcode":<20}: {self.n_bc:}',
            f'{"Unique sample name":<20}: {self.name_is_unique}',
            f'{"Unique sub name":<20}: {self.sub_name_is_unique}',
            f'{"No. of reads":<20}: {self.n_reads:,} M',
            f'{"No. of bases":<20}: {self.n_reads*0.3:,} G (PE150)',
            '='*80,
        ])
        print(msg)
        # incase duplicated names
        if not self.name_is_unique:
            print('> Duplicated sample_name:')
            print(self.db[self.db.duplicated(subset=['name'])])


def get_args():
    """
    extract index table from Excel sample_sheet file
    output:
    id,name,i7_name,i5_name,bc_name,i7,i5,bc,readsM
    YY278s001,CnT_hela_pcf11_sc514158_rep1,Next_Ad2.1,NULL,NULL,TAAGGCGA,NULL,NULL,15
    """
    example = '\n'.join([
        '$ python sample_sheet.py -o yy278.demx.csv -i yy278_sample_sheet.xslx',
    ])
    parser = argparse.ArgumentParser(
        prog='SampleSheet',
        description='Parse sample_sheet',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--excel-file', dest='excel_file', required=True,
        help='sample table in xlsx format, eg: YY00.xlsx')
    parser.add_argument('-o', '--out-csv', dest='out_csv',
        help='save index to csv file')
    parser.add_argument('-f', '--fmt', dest='format', type=int, default=2,
        help='output format, 1=version1, 2=version2, default: [2]')
    return parser


def main():
    args = vars(get_args().parse_args())
    SampleSheet(**args).run()


if __name__ == '__main__':
    main()

#