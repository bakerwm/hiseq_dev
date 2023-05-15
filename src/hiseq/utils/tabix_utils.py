#!/usr/bin/env python3
"""
Utils for tabix 

1. create tabix
2. retrieve data from indexed file (BED, GFF, VCF, SAM, ...)
3. list chromosome from tabix file

modules:

class() Tabix
list_chroms() : tabix --list-chroms 
fetch(region=) : 
fetch(region_list=) :

is_valid_query(query=): 
format_query(query=):


tabix is part of HTSlib
see: https://github.com/samtools/htslib 
python: https://github.com/slowkow/pytabix
file_format: https://github.com/samtools/hts-specs 

Cite:
HTSlib: C library for reading/writing high-throughput sequencing data
James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies
GigaScience, Volume 10, Issue 2, February 2021, giab007, https://doi.org/10.1093/gigascience/giab007

Install tabix

$ wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
$ tar -jxvf htslib-1.17.tar.bz2
$ cd htslib-1.17
$ autoreconf -i 
$ ./configure --prefix ~/.local
$ make
$ make install

Install pytabix
$ pip install pytabix

Create tabix for BED file
$ sort -k1,1 -k2,2n -o demo.bed demo.bed
$ bgzip demo.bed
$ tabix -p bed demo.bed.gz
# create a index file: demo.bed.gz.tbi

List chromosome
$ tabix -l demo.bed.gz 

Fetch records from specific region
$ tabix demo.bed.gz chr1:3205900-4497354
chr1    3205900 3671498 Xkr4    255     -       protein_coding
chr1    3999556 4409241 Rp1     255     -       protein_coding
chr1    4490930 4497354 Sox17   255     -       protein_coding
"""

import os
import sys
import re
import tabix # pytabix
import logging
from subprocess import Popen, PIPE
from itertools import repeat
from multiprocessing import Pool


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


class Tabix(object):
    """    
    Example:
    >>> x = 'demo.bed.gz'
    >>> tbi = Tabix(x)
    """
    def __init__(self, x, **kwargs):
        self.x = x # bgzipped file, or openned bgzipped file
        self.init_args()


    def init_args(self):
        tag = False
        if isinstance(self.x, str):
            if not os.path.exists(self.x):
                log.warning(f'x= file not exists, {self.x}')
            tbi = self.x+'.tbi'
            if self.x.endswith('.gz') and os.path.exists(tbi):
                tag = True
            else:
                x_name = os.path.splitext(self.x)[0]
                log.warning(f'tabix file not detected, {tbi}')
                log.warning(
                '''
                # Prepare a tabix file with the following commands
                # 1. sort the BED file
                sort -k1,1 -k2,2n -o demo.bed demo.bed
                # 2. bgzip compress file
                bgzip demo.bed
                # 3. Create tabix
                tabix -p bed demo.bed.gz
                '''
                )
        else:
            log.error(f'invalid x=, expect str, got {type(self.x).__name__}')
        self.fh = self.open(self.x)


    def open(self, x):
        """
        Open a bgzipped and tabix indexed file with pytabix
        alternative: shell command tabix
        use subprocess

        Parameters:
        x (str): path to a bgzipped file, including tabix in the same folder

        Returns:
        tabix.open object
        """
        try:
            out = tabix.open(x)
        except:
            out = None
            log.error(f'Could not open gzipped file: {x}')
        return out


    def format_region(self, region=None, chr=None, start=None, end=None):        
        if self.is_valid_region(region):
            out = region
        elif isinstance(chr, str):
            if isinstance(start, int) and isinstance(end, int):
                start, end = [abs(start), abs(end)]
                if start > end:
                    start, end = [end, start] # make sure start < end
                out = f'{chr}:{start}-{end}'
            else:
                out = chr
        else:
            out = None
        return out if self.is_valid_region(out) else None


    def is_valid_region(self, region):
        """
        Format query string as: chr:start-end for tabix 
        1. chr1:1-100
        2. chr1
        """
        p1 = re.compile('^(\w+):(\d+)-(\d+)$|^(\w+)$')
        if isinstance(region, str):
            out = isinstance(p1.match(region), re.Match)
        else:
            out = False
        return out


    def list_chroms(self):
        """
        List all chromosome names
        see:
        tabix -l demo.bed.gz
        """
        p = Popen(['tabix', '-l', self.x], stdout=PIPE)
        return [i.decode().strip() for i in p.stdout]


    def fetch(self, region=None, chr=None, start=None, end=None):
        """
        Fetch records by chr:start-end from tabix indexed file
        
        Parameters:
        fh (tabix.open) a tabix.open object
        pos (str): position to extract, format: chr:start-end 
        chr (str): chromosome name 
        start (int): start 
        end (int): end

        Examples:
        tabix.open(x).query("chr1", 1, 1000)
        tabix.open(x).queryi(0, 1, 1000)
        tabix.open(x).querys("chr1:1-1000")
        see: https://github.com/slowkow/pytabix

        Returns:
        records
        """
        if isinstance(self.fh, tabix.open):
            query = self.format_region(region, chr, start, end) # chr1:1-100, chr1
            if isinstance(query, str):
                try:
                    out = self.fh.querys(query) # tabix.iter
                except:
                    out = None # show error ?
                return out


    def fetch2(self, region_list):
        for i in region_list:
            if self.is_valid_region(i):
                yield self.fetch(i) # generator


    def fetch_by_bed(self, bed):
        if isinstance(bed, str):
            if os.path.exists(bed):
                with open(bed) as r:
                    for l in r:
                        s = l.strip().split() # tab
                        region = f's[0]:s[1]-s[2]'
                        yield self.fetch(region)


