#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert sample_sheet.xlsx to csv files
sheet_name: sample_sheet
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
import pathlib
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
        self.df1 = self.load_hiseq_index() # name:seq
        self.df2 = {v:k for k,v in self.df1.items()} # seq:name
        if isinstance(x, str):
            self.name = self.get_name(x)
            self.index = self.get_index(x)
            self.is_valid = self.is_name(x) or self.is_index(x)
            # x in self.df1or x in self.df2
        elif isinstance(x, list):
            self.name = [self.get_name(i) for i in x]
            self.index = [self.get_index(i) for i in x]
            self.is_valid = [self.is_name(i) or self.is_index(i) for i in x]
        else:
            raise ValueError('illegal x={}, expect str,list got {}'.format(
                x, type(x).__name__))


    def load_hiseq_index(self):
        with res.path('hiseq.data', 'illumina_index.csv') as f:
            ff = str(f)
        d = {}
        with open(ff) as r:
            for line in r:
                if line.startswith('#'):
                    continue
                name,seq = line.strip().split(',')
                d.update({name:seq})
        return d


    def is_name(self, x):
        return x in self.df1 or x.upper() == 'NULL'


    def is_index(self, x):
        return x in self.df2 or x.upper() == 'NULL'


    def get_name(self, x):
        if self.is_name(x):
            out = x
        elif self.is_index(x):
            out = self.df2.get(x, 'NULL')
        else:
            out = 'NULL'
        return out


    def get_index(self, x):
        if self.is_index(x):
            out = x
        elif self.is_name(x):
            out = self.df1.get(x, 'NULL')
        else:
            out = 'NULL'
        return out


class SampleSheet(object):
    """
    Processing sample table file.xlsx
    ## Example
    >>> s = SampleSheet(x='YY00.xlsx', out_dir='data')
    >>> s.to_MGI_table()
    >>> s.to_barcode_table()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        if self.x.endswith('.xlsx'):
            self.df = self.load_xlsx(self.x)
        elif self.x.endswith('.csv'):
            self.df = self.load_csv(self.x)
        else:
            raise ValueError('illegal x={}'.format(self.x))


    def init_args(self):
        args_init = {
            'x': None,
            'out_dir': None,
        }
        self = update_obj(self, args_init, force=False)
        if not os.path.isfile(self.x):
            raise ValueError('file not exists: {}'.format(self.x))
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        # try...except... 
        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir)
        self.out_dir = os.path.abspath(self.out_dir)
        self.x = os.path.abspath(self.x)
        self.xlsx_prefix = os.path.splitext(os.path.basename(self.x))[0]
        self.mgi_csv = os.path.join(self.out_dir, self.xlsx_prefix + '.MGI.csv')
        self.demx_csv = os.path.join(self.out_dir, self.xlsx_prefix + '.demx.csv')


    def load_csv(self, x):
        """
        Sample sheet in csv format:
        name, i7_id, i5_id, bc_id, reads
        output: (index_name)
        """
        c = ['name', 'i7', 'i5', 'bc', 'reads']
        try:
            df = pd.read_csv(x, header=None, comment='#')
            df.fillna('NULL', inplace=True) # convert NaN to 'NULL'
            n_rows, n_cols = df.shape
            if n_cols == 4: #missing reads
                df.columns = c[:4]
                df = df.assign(reads=0)
            elif n_cols == 5:
                df.columns = c
            else:
                raise ValueError('expect 4|5 columns in csv, got {}'.format(n_cols))
            # convert to index_name
            df = df.assign(
                i7 = HiSeqIndex(df['i7'].to_list()).name,
                i5 = HiSeqIndex(df['i5'].to_list()).name,
                bc = HiSeqIndex(df['bc'].to_list()).name,
            )
        except ValueError as e:
            print(e)
            df = pd.DataFrame(columns=c)
        return df # ['name', 'i7', 'i5', 'bc', 'reads']


    def load_xlsx(self, x):
        c = ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']
        try:
            df = pd.read_excel(x, sheet_name='sample_sheet')
            df = df.dropna(thresh=2)
            df = df[c]
            df.columns = ['name', 'i7', 'bc', 'reads']
            df['name'] = self.sanitize(df['name'].to_list()) # sanitize
            df['i5'] = 'NULL'
        except:
            df = pd.DataFrame(columns=['name', 'i7', 'i5', 'bc', 'reads'])
        return df # ['name', 'i7', 'i5', 'bc', 'reads']


    def sanitize(self, x):
        """
        The filename rules:
        1. [A-Za-z0-9_-.]
        2. convert to '_'
        3. unique
        """
        if isinstance(x, str):
            out = re.sub('[^A-Za-z0-9_\-]', '_', x)
            out = re.sub('_+', '_', out)
        elif isinstance(x, list):
            out = [self.sanitize(i) for i in x]
        else:
            log.warning('illegal x')
            out = x
        return out


    def to_demx_table(self):
        """
        Convert to demx table
        format:
        sample_name,i7,i5,barcode
        """
        df = self.df # name, i7, bc, reads
        # convert to index_seq
        df = df.assign(
            i7 = HiSeqIndex(df['i7'].to_list()).index,
            i5 = HiSeqIndex(df['i5'].to_list()).index,
            bc = HiSeqIndex(df['bc'].to_list()).index,
        )
        df = df.loc[:, ['name', 'i7', 'i5', 'bc', 'reads']]
        df.to_csv(self.demx_csv, index=False, header=False)
        return df


    def to_MGI_table(self):
        """
        Convert to MGI table
        1. Assign the names by i7_index_id
        2. concatenate files with same P7 index
        3. Sum the reads
        """
        df2 = self.df.groupby('i7').sum() # required: i7, reads
        gb = df2['reads']*0.3
        df2 = df2.assign(
            i7_seq = HiSeqIndex(df2.index.to_list()).index,
            GB = gb.round(1),
        )
        df2.reset_index(inplace=True)
        # add barcode width
        df2 = df2.assign(
            index = self.sanitize(df2['i7'].to_list()),
            i7_width = df2['i7_seq'].apply(len),            
        )
        df2 = df2.loc[:, ['index', 'i7_seq', 'i7_width', 'GB', 'reads']]
        df2.to_csv(self.mgi_csv, index=False)
        return df2


    def to_barcode_table(self):
        """
        Generate index table for barcode:
        format:
        filename,i7,i5,bc
        sample1,NULL,NULL,CCTATA
        filename:
        i7_index.csv
        """
        bc_tables = []
        s = self.df.groupby('i7').size()
        s = s[s>1]
        i7_list = s.index.to_list()
        if len(i7_list) > 0:
            for i7 in i7_list:
                i7_csv = os.path.join(self.out_dir, 'i7_index.{}.csv'.format(i7))
                df1 = self.df.loc[self.df['i7'] == i7]
                # convert name to index
                df2 = df1.assign(
                    bc_seq = HiSeqIndex(df1['bc'].to_list()).index,
                    i7_seq = 'NULL', # HiSeqIndex(df1['i7'].to_list()).index,
                    i5_seq = 'NULL',
                )
                df2 = df2.loc[:, ['name', 'i7_seq', 'i5_seq', 'bc_seq']]
                df2.to_csv(i7_csv, index=False, header=False)
                bc_tables.append(i7_csv)
        return bc_tables


    def run(self):
        s = self.df.groupby('i7').size()
        s = s[s>1]
        # message
        df_demx = self.to_demx_table()
        df_mgi = self.to_MGI_table()
        bc_tables = self.to_barcode_table()
        n_smp = self.df.shape[0]
        n_i7 = df_mgi.shape[0]
        n_i7_bc = len(s)
        n_total = self.df['reads'].sum()
        msg = '\n'.join([
            '='*80,
            '{:<20}: {:}'.format('Program', 'SampleSheet'),
            '{:<20}: {:}'.format('Date', get_date()),
            '{:<20}: {:}'.format('Input table', self.x),
            '{:<20}: {:}'.format('to MGI table', self.mgi_csv),
            '{:<20}: {:}'.format('to Demx table', self.demx_csv),
            '{:<20}: {:}'.format('No. of samples', n_smp),
            '{:<20}: {:}'.format('No. of i7', n_i7),
            '{:<20}: {:}'.format('No. of i7 with bc', n_i7_bc),
            '{:<20}: {:,}M'.format('No. of reads', int(n_total)),
            '{:<20}: {:,}G'.format('No. of bases', int(n_total*0.3)),
            '='*80,
        ])
        print(msg)


def get_args():
    """
    Prepare sample sheet for Demx
    output:
    1. sample_name,i7,i5,barcode
    2. i7_name,i7,reads
    3. sample_name,NULL,NULL,barcode (bc only)
    """
    example = '\n'.join([
        'input sample_sheet.xslx',
        '1. parse sample_sheet',
        '$ python sample_sheet.py -s YY00.xslx -o data',
    ])
    parser = argparse.ArgumentParser(
        prog='SampleSheet',
        description='Parse sample_sheet',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--xlsx-table', dest='x', required=True,
        help='sample table in xlsx format, eg: YY00.xlsx')
    parser.add_argument('-o', '--out-dir', dest='out_dir',
        help='directory to save the results')
    return parser
        
        
def main():
    args = vars(get_args().parse_args())
    SampleSheet(**args).run()


if __name__ == '__main__':
    main()

#