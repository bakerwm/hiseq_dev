#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Help functions for hiseq module
# function
is_supported
is_hiseq_dir
list_hiseq_file
list_hiseq_dir
read_hiseq
# class
HiseqReader
"""


import os
# import importlib
import shutil
import glob
import importlib.resources as res
from hiseq.utils.utils import log, update_obj, Config #, run_shell_cmd
# from hiseq.utils.file import file_exists, check_dir
from hiseq.utils.seq import list_fx, fx_name, check_fx_paired
from hiseq.utils.file import list_dir, file_exists


## functions for hiseq-global
def is_supported(x=True, key=True, return_values=False):
    """
    Check if the aligner/genome is supported or not
    Parameters
    ----------
    x : str, bool
        if True, match all values
        
    key : str, bool
        options: ['aligner', 'genome']
        default: [True]
        
    return_values :  bool
        Return the specific value, by 'key', ignore 'x'
    
    Check args(x) is supported: 
      - aligner
      - genomes      
    Config file saved in: hiseq/data/config.yaml
    Examples:
    1. list all supported aligners
    >>> is_supported(key='aligner', return_values=True) # [...]
    2. list all supported genomes
    >>> is_supported(key='genome', return_values=True) # [...]
    3. check if specific aligner/genome is supported
    >>> is_supported('dm6') # True
    >>> is_supported('dm6', 'aligner') # False
    >>> is_supported('bowtie2', 'aligner') # True
    >>> is_supported('bowtie2', 'genome') # False
    >>> is_supported('abc') # False
    # check first:
    1. key is True or valid str.
    2. x is valid or not
    3. return all values
    # how-to
    1. x is True. (keys)
    2. x is True, key is True. (all key:values) 
    3. x is str, key is True. (...) 
    4. x is str, key is str. (...)
    """
    # config file saved
    with res.path('hiseq.data', 'config.yaml') as f:
        cfg = str(f)
    d = Config().load(cfg)
    out = False
    if not isinstance(d, dict):
        log.error('could not read config: {}'.format(cfg))
        return False
    if x is True:
        if key is True:
            out = d if return_values else True
        elif key in d:
            out = d.get(key, False) if return_values else True
        else:
            out = False
    elif isinstance(x, str):
        if key is True:
            out = x in [j for i in list(d.values()) for j in i]
        elif key in d:
            out = x in d.get(key, False)
        else:
            out = False
        if return_values and out:
            out = x
    else:
        out = False
    return out


def is_hiseq_dir(x, hiseq_type='auto'):
    """
    Check if the input is hiseq dir or not
    x could be str, list
    """
    if isinstance(x, str):
        a = read_hiseq(x)
        k = False # init
        if a.is_hiseq:
            if hiseq_type == 'auto' or hiseq_type is True:
                k = True
            elif isinstance(hiseq_type, str):
                k = any([
                    a.hiseq_type == hiseq_type,
                    a.hiseq_type.startswith(hiseq_type),
                    a.hiseq_type.endswith(hiseq_type),
                ])
    elif isinstance(x, list):
        k = [is_hiseq_dir(i, hiseq_type) for i in x]
    else:
        log.error('x is {}, expect str, list'.format(type(x).__name__))
        k = None
    return k


def list_hiseq_file(x, keys='bam', hiseq_type='auto'):
    """
    Return the specific file from project_dir of hiseq
    Parameters
    ----------
    x:  str
        Directory of the hiseq project
    keys:  str
        keys of the file in hiseq project, eg: 'bam', 'align_dir', ...
        default: ['bam']
    hiseq_type:  str
        Type of the hiseq, ['r1', 'rn', 'rx', 'rt', 'rp']
        see: HiseqReader() for details; default: ['r1']
    """
    out = []
    dirs = list_hiseq_dir(x, hiseq_type) # all dirs
    if dirs is not None:
        for i in dirs:
            a = read_hiseq(i)
            out.append(getattr(a, keys, None))
    # to str !!!
    if len(out) == 1:
        out = out[0]
    return out


def list_hiseq_dir(x, hiseq_type='auto'):
    """
    Return the project_dir of hiseq
    Parameters
    ----------
    x:  str
        Directory of the hiseq project
    hiseq_type:  str
        Type of the hiseq, ['r1', 'rn', 'rx', 'rt', 'rp']
        see: HiseqReader() for details; default: ['auto'],
    """
    a = read_hiseq(x)
    out = []
    if a.is_hiseq:
        # update hiseq_type
        if isinstance(hiseq_type, str):
            if(hiseq_type == 'auto'):
                hiseq_type = a.hiseq_type # auto
        if any([a.hiseq_type.endswith(i) for i in ['r1', 'rp', 'merge', 'alignment']]):
            out = [x]
        elif a.hiseq_type.endswith('rn'):
            out = [x] # add rn
            r1 = getattr(a, 'rep_list', None)
            if isinstance(r1, list):
                out += r1 # add r1
        elif a.hiseq_type.endswith('rx'): 
            if a.hiseq_type.startswith('atac_'):
                out = list_dir(x, include_dir=True)
                # out = [i for i in out if is_hiseq_dir(i)]
            else:
                # for rn dirs
                rn = []
                if any([a.hiseq_type.startswith(i) for i in ['cnr', 'cnt', 'chip']]):
                    kt = ['ip_dir', 'input_dir']
                elif a.hiseq_type.startswith('rnaseq_'):
                    kt = ['wt_dir', 'mut_dir']
                else:
                    kt = []
                for k in kt:
                    kd = getattr(a, k, None)
                    if isinstance(kd, str):
                        rn += [kd]
                # for r1
                try:
                    r1 = [j for i in rn for j in list_hiseq_dir(i, 'r1')]
                except:
                    r1 = []
                    # print('!B-1', rn)
                out = rn + r1
                # for rx
                out.append(x)
                # out = [i for i in out if is_hiseq_dir(i)]
        else:
            out = [] # empty        
        out = list(set(out)) # remove duplicates
        out = [i for i in out if is_hiseq_dir(i, hiseq_type)] # keep hiseq
        # fix None
        out = [i for i in out if i is not None]
#         if out is None:
#             out = []
    return sorted(out)


def read_hiseq(x, hiseq_type='auto'):
    """
    Parameters
    ---------
    x:  str
        Path to Hiseq directory
    hiseq_type: str
        Check the x is one-of-hiseq_types, could be head/tail;
        atacseq_r1, r1, rn, rnaseq_rx, ...
        default: [True], do
    Read config from hiseq directory
    """
    a = HiseqReader(x)
    if hasattr(a, 'hiseq_type'):
        a_hiseq_type = a.hiseq_type
        k = False
        if hiseq_type is True:
            k = True #pass
        elif isinstance(hiseq_type, str):
            k = any([
                a_hiseq_type == hiseq_type,
                a_hiseq_type.startswith(hiseq_type),
                a_hiseq_type.endswith(hiseq_type),
            ])
        else:
            pass
#         # check
#         if not k:
#             log.warning('hiseq_dir not match, expect {}, got {}'.format(
#                 hiseq_type, a_hiseq_type))
    else:
        # raise ValueError('not a hiseq dir: {}'.format(x))
        # log.warning('not a hiseq dir: {}'.format(x))
        pass # reduce messagage #
    return a


class HiSeqReader(object):
    """
    Parameters
    ---------
    x: str
        Path to Hiseq directory
    """
    def __init__(self, x):
        self.x = x
        self.init_args()


    def init_args(self):
        self.is_hiseq = False # init
        d = self.load()
        if isinstance(d, dict):
            hiseq_types = [
                'atacseq_type', 'rnaseq_type', 'hiseq_type',
                'align_type'
            ] # !!!! to-be-update
            # which hiseq
            hiseq_type = None
            for h in hiseq_types:
                if h in d:
                    hiseq_type = h
                    break
            # which type
            if isinstance(hiseq_type, str):
                a = d.get(hiseq_type, None)
            else:
                a = None
            # check output
            if isinstance(a, str):
                self.is_hiseq = True
                self.hiseq_type = a
                self.is_hiseq_r1 = a.endswith('_r1') # fq
                self.is_hiseq_rn = a.endswith('_rn') # group
                self.is_hiseq_rx = a.endswith('_rx') # group vs group
                self.is_hiseq_rt = a.endswith('_rt') # !!
                self.is_hiseq_rp = a.endswith('_rp') # report
                self.is_hiseq_merge = a.endswith('_merge') # merge multiple dirs
#             else:
#                 log.warning('unknown hiseq dir: {}'.format(self.x))
        else:
            # log.warning('not a hiseq dir: {}'.format(self.x))
            pass # reduce message log
        # update args
        self.args = update_obj(self, d, force=True)


    def list_config(self):
        """
        List the config files
        Support: yaml, pickle, toml, json, ... [priority]
        return file list
        # hiseq
        hiseq
          |-config
          |   |-config.yaml

        # alignment
        align_dir
          |- smp_nmae
          |    |- index
          |    |    |- config.pickle
        """
        c_files = ['config.' + i for i in ['yaml', 'pickle', 'toml', 'json']]
        # search config files
        out = None
        for f in c_files:
            c1 = os.path.join(self.x, 'config', f)
            c2 = os.path.join(self.x, '*', '*', f)
            c1x = glob.glob(c1)
            c2x = glob.glob(c2)
            if len(c1x) > 0:
                out = c1x[0]
                break
            elif len(c2x) > 0:
                out = c2x[0]
                break
            else:
                continue
        return out


    def load(self):
        if isinstance(self.x, str) and os.path.exists(self.x):
            config_file = self.list_config()
            out = Config().load(config_file)
        else:
            out = None
        return out


class HiSeqDesignAtac(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)


    def init_args(self):
        args_init = {
            'design': None,
            'fq_dir': None,
            'fq1': None,
            'fq2': None,
            'append': True,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'atac_rd'
        self.build_design = False # force


    def list_fq_files(self, x):
        if isinstance(x, str):
            f_list = list_fx(x)
        if len(f_list) < 2:
            sub_dirs = list_dir(x, include_dirs=True)
            sub_files = [list_fx(i) for i in sub_dirs if os.path.isdir(i)]
            f_list = [f for i in sub_files if isinstance(i, list) for f in i]
        if len(f_list) < 2:
            log.error('not enough fq files: {}'.format(x))
            return None
        return f_list

    
    def parse_fq(self):
        if isinstance(self.fq_dir, str):
            f_list = self.list_fq_files(self.fq_dir)
        elif isinstance(self.fq1, list) and isinstance(self.fq2, list):
            f_list = self.fq1 + self.fq2
        elif isinstance(self.fq1, str) and isinstance(self.fq2, str):
            f_list = [self.fq1, self.fq2]
        else:
            f_list = []
        if len(f_list) < 2: # paired
            log.error('no fq files found, check -r, -1, -2')
            return None
            # raise ValueError('fq_dir, fq1,fq2, failed')
        # separate fq1, fq2
        fq1 = [i for i in f_list if fx_name(i).endswith('1')]
        fq2 = [i for i in f_list if fx_name(i).endswith('2')]
        if not check_fx_paired(fq1, fq2):
            log.error('fq not paired, check -r, -1, -2')
            return None
        # recognize samples by name: filename_repx.fq.gz
        d = {}
        f_names = fx_name(f_list, fix_pe=True, fix_rep=True, fix_unmap=True)
        f_names = sorted(list(set(f_names)))
        for f in f_names:
            df = d.get(f, {})
            for i,j in zip(fq1, fq2):
                iname1 = fx_name(i, fix_pe=True, fix_rep=True, fix_unmap=True)
                if iname1 == f:
                    iname2 = fx_name(i, fix_pe=True, fix_rep=False, fix_unmap=True)
                    df.update({
                        iname2: [i, j]
                    })
            # update
            d.update({f: df})
        return d


    def update_design(self):
        # load design/fq_groups
        if file_exists(self.design) and self.append:
            d = Config().load(self.design)
        else:
            d = {}
        # new fq_groups
        d.update(self.parse_fq())
        return d


    def show_msg(self):
        msg = ['='*80]
        n_fq = 0
        for k,v in self.fq_groups.items():
            # smp_name, ...
            msg.append('>{}'.format(k))
            for i,j in v.items():
                n_fq += 1
                # rep1, ...
                msg.append(
                    '\n'.join([
                        '{} {}'.format(' '*2, i),
                        '{} r1: {}'.format(' '*4, j[0]),
                        '{} r2: {}'.format(' '*4, j[1]),
                    ])
                )
        msg.append('\n'.join([
            '{:>8}: {}'.format('samples', len(self.fq_groups)),
            '{:>8}: {}'.format('fq pairs', n_fq)
        ]))
        msg.append('='*80)
        out = '\n'.join(msg)
        print(out)


    def run(self):
        self.fq_groups = self.parse_fq()
        self.show_msg()
        Config().dump(self.fq_groups, self.design)


# def main():
#     args = {
#         'design': 'aaa.yaml',
#         'fq_dir': '/data/yulab/wangming/work/wmlib/hiseq_dev/test/data',
#     }
#     HiSeqDesignAtac(**args).run()


# if __name__ == '__main__':
#     main()
