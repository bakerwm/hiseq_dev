#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genome related files
files:
  - fasta 
  - fa_size
  - index [bowtie, bowtie2, star, bwa, ...]
  - bed 
  - gtf
  - te (gtf)
  - piRC (gtf)
  - blacklist  
check:
  - supported (hiseq.utils.is_supported)
"""


import os
import shutil
import pathlib
import argparse
import yaml
import pysam
import importlib.resources as res
from xopen import xopen
from hiseq.utils.utils import log, update_obj, Config, download_file
from hiseq.utils.file import file_exists, file_abspath, check_dir
from hiseq.utils.hiseq_utils import is_supported


class SetupGenome(object):
    """
    retrieve genome data from dir
    update genome data from file (yaml)
    template
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        gd = os.path.join(str(pathlib.Path.home()), 'data', 'genome')
        self.genome_dir = getattr(self, 'genome_dir', None)
        if self.genome_dir is None:
            self.genome_dir = gd
        self.genome_dir = file_abspath(self.genome_dir)
        self.overwrite = getattr(self, 'overwrite', False)


    def template(self, out=None):
        d = {
            'genome': {
                'fasta': None,
                'fai': None,
                'gtf': {
                    'ensembl': None,
                    'refseq': None,
                    'genecode': None,
                },
                'bed': {
                    'ensembl': None,
                    'refseq': None,
                    'gencode': None,
                },
                'te': None,
                'piRC': None,
                'blacklist': None,
                'phylop100': None,
            }
        }
        if isinstance(out, str):
            out_dir = os.path.dirname(os.path.abspath(out))
            check_dir(out_dir) # create out_dir
            if os.path.exists(out):
                log.warning('file exists: {}'.foramt(out))
            else:
                Config().dump(d, out)
        else:
            # to STDOUT
            print(yaml.dump(d, default_flow_style=False))
        return d


    def get_fa(self, genome):
        """
        genome:
        genome_dir: 
        """
        if is_supported(genome, key='genome'):
            fa1 = os.path.join(self.genome_dir, genome, 'bigZips', genome+'.fa')
            fa2 = os.path.join(self.genome_dir, genome, 'fasta', genome+'.fa')
            out = fa1 if file_exists(fa1) else fa2 if file_exists(fa2) else None
        else:
            out = None
        return out


    def get_fai(self, genome):
        fa = self.get_fa(genome)
        out = None
        if isinstance(fa, str):
            fai1 = fa + '.fai'
            fai2 = os.path.splitext(fa)[0] + '.chrom.sizes'
            if file_exists(fai2):
                out = fai2
            elif file_exists(fai1):
                out = fai1
            else:
                pysam.faidx(fa)
                out = fai1
        return out


    def get_gene_anno(self, genome, fmt='bed', db='ensembl', rmsk=False):
        """
        Return the gene annotation in BED,GTF format
        support UCSC, ensembl, gencode
        """
        if is_supported(genome, key='genome'):
            ad = os.path.join(self.genome_dir, genome, 'annotation_and_repeats')
            anno_dir = getattr(self, 'anno_dir', ad)
            if rmsk:
                ext = '.{}.rmsk.{}'.format(db, fmt)
            else:
                ext = '.{}.{}'.format(db, fmt)
            anno = os.path.join(anno_dir, genome + ext)
            out = anno if file_exists(anno) else None
        else:
            out = None
        return out


    def get_te(self, genome, fmt='gtf'):
        """
        Return TE annotation of the genome (dm6)
        or return TE consensus sequence for the genome (dm3)
        """
        if is_supported(genome, key='genome'):
            tn = '{}_transposon.{}'.format(genome, fmt)
            td = os.path.join(self.genome_dir, genome, genome + '_transposon')
            te = os.path.join(td, tn)
            out = te if file_exists(te) else None
        else:
            out = None
        return out

    
    def get_piRC(self, genome, fmt='gtf'):
        """
        Return piRNA cluster annotation of the genome (dm6)
        or return piRNA cluster consensus sequence for the genome (dm6)
        """
        if is_supported(genome, key='genome'):
            pn = '{}_piRNA_clusters.{}'.format(genome, fmt)
            pd = os.path.join(self.genome_dir, genome, genome + '_piRNA_clusters')
            pc = os.path.join(pd, pn)
            out = pc if file_exists(pc) else None
        else:
            out = None
        return out


    def get_phylop100(self, genome):
        """
        Return the phylop100 bigWig file of hg19, only
        for conservation analysis
        """
        if is_supported(genome, key='genome'):
            pp = os.path.join(
                self.genome_dir, genome, 'phyloP100way',
                genome + '.100way.phyloP100way.bw'
            )
            out = pp if file_exists(pp) else None
        else:
            out = None
        return out


    def get_blacklist(self, genome):
        """
        blacklist files were downloaded from github repo:
        https://github.com/Boyle-Lab/Blacklist/
        supported:
        dm3, dm6, ce10, ce11, mm10, hg19, and hg38        
        the local path is:
        genome_dir/{genome}/annotation_and_repeats/blacklist/{genome}.blacklist.v2.bed
        """
        if not is_supported(genome, key='genome'):
            return None
        # remote files
        url = 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists'
        fn = {
            'dm3': 'dm3-blacklist.v2.bed.gz',
            'dm6': 'dm6-blacklist.v2.bed.gz',
            'ce10': 'ce10-blacklist.v2.bed.gz',
            'cd11': 'cd11-blacklist.v2.bed.gz',
            'mm10': 'mm10-blacklist.v2.bed.gz',
            'hg19': 'hg19-blacklist.v2.bed.gz',
            'hg38': 'hg38-blacklist.v2.bed.gz'
        }
        fd = {k:os.path.join(url, v) for k,v in fn.items()} # remote
        # update path
        if not genome in fn:
            return None
        # local file
        ad = os.path.join(self.genome_dir, genome, 'annotation_and_repeats')
        anno_dir = getattr(self, 'anno_dir', ad)
        fc = {
            k:os.path.join(
                anno_dir, 'blacklist', os.path.splitext(v)[0]
            ) for k,v in fn.items()
        }
        out = fc.get(genome, None)
        if not file_exists(out):
            # download from github
            f_src = fd.get(genome, None)
            try:
                out_gz = out + '.gz' # gzipped
                check_dir(os.path.dirname(out))
                download_file(f_src, out_gz) # gzipped
                with xopen(out_gz, 'rb') as r, xopen(out, 'wb') as w:
                    shutil.copyfileobj(r, w)
            except:
                msg = '\n'.join([
                    'failed downloading blacklist',
                    'Using the following command in terminal to download file manually:',
                    '$ mkdir -p {}'.format(os.path.dirname(out)),
                    '$ wget -O {} {}'.format(out+'.gz', f_src),
                    '$ gunzip {}'.format(out+'.gz')
                ])
                log.error(msg)
                out = None
        return out


    def is_valid(self, d):
        """
        required: fasta, gtf, bed
        """
        fa = d.get('fasta', None)
        f1 = file_exists(fa)
        gtf = d.get('gtf', None) # db
        if isinstance(gtf, dict):
            f2 = any([file_exists(i) for i in list(gtf.values())])
        else:
            f2 = False
        bed = d.get('bed', None) # db
        if isinstance(bed, dict):
            f3 = any([file_exists(i) for i in list(bed.values())])
        else:
            f3 = False
        return all([f1, f2, f3])


    def run(self):
        if not file_exists(self.genome_dir):
            log.error('genome_dir not exists: {}'.format(self.genome_dir))
            return None
        # global
        with res.path('hiseq.data','genome.yaml') as f: # importlib.resources
            gf = str(f)
        df = Config().load(gf) # default
        if df is None:
            df = {} # init
        # check available genome
        n_update = 0
        for g in is_supported(key='genome', return_values=True):
            # latest, local values
            val = {
                'fasta': self.get_fa(g),
                'fai': self.get_fai(g),
                'gtf': {
                    'ensembl': self.get_gene_anno(g, fmt='gtf', db='ensembl'),
                    'refseq': self.get_gene_anno(g, fmt='gtf', db='refseq'),
                    'gencode': self.get_gene_anno(g, fmt='gtf', db='gencode'),
                },
                'bed': {
                    'ensembl': self.get_gene_anno(g, fmt='bed', db='ensembl'),
                    'refseq': self.get_gene_anno(g, fmt='bed', db='refseq'),
                    'gencode': self.get_gene_anno(g, fmt='bed', db='gencode'),
                },
                'te': self.get_te(g, 'gtf'),
                'piRC': self.get_piRC(g, 'gtf'), 
                # 'blacklist': self.get_blacklist(g),
                'phylop100': self.get_phylop100(g),
            }
            if self.is_valid(val):
                val.update({'blacklist': self.get_blacklist(g)})
                if g in df and not self.overwrite:
                    pass
                else:
                    df.update({g:val})
                n_update += 1
        # filt all config
        dg = {k:v for k,v in df.items() if self.is_valid(v)}
        # report msg
        msg = '\n'.join([
            '#'*80,
            '{:>10} : {}'.format('program', 'Setup Genome'),
            '{:>10} : {}'.format('genome_dir', self.genome_dir),
            '{:>10} : {}'.format('update', n_update),
            '{:>10} : {}'.format('genome', ','.join(list(dg.keys()))),
            '{:>10} : {}'.format('config', gf),
            '{:>10} : {}'.format('status', 'pass' if len(dg) > 0 else 'failed'),
            '#'*80,
        ])
        print(msg)
        # print(yaml.dump(d, default_flow_style=False))
        if len(dg) > 0:
            Config().dump(dg, gf)
        else:
            log.error('No genome found in: {}'.format(self.genome_dir))


class Genome(object):
    """
    Retrieve genome information:
    fasta, fai, gtf, bed, te, piRC, blacklist
    """
    def __init__(self, genome, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.genome = genome
        db = self.load_db() # all
        gv = db.get(genome, {})
        # db_src: ensembl, refseq, gencode
        db_src = getattr(self, 'db_src', 'ensembl')
        fmt = getattr(self, 'fmt', 'gtf')
        self.fasta = gv.get('fasta', None)
        self.fai = gv.get('fai', None)
        self.te = gv.get('te', fmt)
        self.piRC = gv.get('piRC', fmt)
        self.blacklist = gv.get('blacklist', None)
        # gtf and bed
        self.gtf = gv.get('gtf', {}).get(db_src, None)
        self.bed = gv.get('bed', {}).get(db_src, None)


    def load_db(self):
        with res.path('hiseq.data', 'genome.yaml') as f:
            gf = str(f)
        df = Config().load(gf)
        if not isinstance(df, dict):
            df = {} # init
        return {k:v for k,v in df.items() if self.is_valid(v)}


    def is_valid(self, d):
        """
        required: fasta, gtf, bed
        """
        if not isinstance(d, dict):
            return False
        fa = d.get('fasta', None)
        f1 = file_exists(fa)
        gtf = d.get('gtf', None) # db
        if isinstance(gtf, dict):
            f2 = any([file_exists(i) for i in list(gtf.values())])
        else:
            f2 = False
        bed = d.get('bed', None) # db
        if isinstance(bed, dict):
            f3 = any([file_exists(i) for i in list(bed.values())])
        else:
            f3 = False
        return all([f1, f2, f3])


def get_args():
    example = '\n'.join([
        'Setup genome:',
        '1. install genomes',
        '$ python genome.py -i genome_dir',
        '2. install genomes, overlap',
        '$ python genome.py -i genome_dir -O',
    ])
    parser = argparse.ArgumentParser(
        prog='setup_genome',
        description='Setup genomes for hiseq',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--genome-dir', dest='genome_dir', required=False,
        help='The dir of genome files')
    parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
        help='overwrite exists genomes')
    return parser


def main():
    args = vars(get_args().parse_args())
    SetupGenome(**args).run()


if __name__ == '__main__':
    main()

# EOF