#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Run fastqc for fastq file, for base-content
set -nogroup (fastqc), --nogroup (falco) for each base
1. only fastq (convert fa to fq)
# first N-base
# last N-base
# reverse-complement
Date: 2022-08-31
"""


import os
import pathlib
import re
import shutil
import argparse
import json
import pandas as pd
from shutil import which
import plotly.express as px
from hiseq.qc.read_fastqc import ReadFastQC
from multiprocessing import Pool
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd
from hiseq.utils.file import fix_out_dir
from hiseq.utils.seq import fx_name


class BaseContentR1(object):
    """
    Check Per Base Content for single fastq file
    using: FastQC or Falco
    """
    def __init__(self, fq, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.fq = fq
        self.init_args()


    def init_args(self):
        args = {
            'out_dir': None,
            'parallel_jobs': 1,
            'overwrite': False,
            'left_end': None,
            'rm_tmp': False,
        }
        self = update_obj(self, args, force=False)
        # if not isinstance(self.out_dir, str):
        #     self.out_dir = str(pathlib.Path.cwd())
        # if not os.path.exists(self.out_dir):
        #     os.makedirs(self.out_dir)
        self.out_dir = fix_out_dir(self.out_dir)
        # name = os.path.basename(self.fq)
        # if name.endswith('.gz'):
        #     name = name[:-3] # remove ".gz"
        # name = os.path.splitext(name)[0] # remove '.fq'
        # self.name = name
        self.name = fx_name(self.fq, fix_pe=True)
        self.out_txt = os.path.join(self.out_dir, self.name+'_fastqc_data.txt')
        self.out_json = os.path.join(self.out_dir, self.name+'_fastqc_data.json')
        self.out_png = os.path.join(self.out_dir, self.name+'_basecontent.png')
        self.plot_json = os.path.join(self.out_dir, self.name+'_basecontent.json')


    def get_df(self):
        if os.path.exists(self.out_txt):
            # for Per Base Content
            m = 'Per base sequence content'
            qc = ReadFastQC(self.out_txt)
            qc.save_as(self.out_json)
            d1 = qc.get_module(m)
            d2 = d1.get('table')
            d2 = {k:list(map(eval, v)) for k,v in d2.items()}
            df = pd.DataFrame(data=d2)
            # subset
            if isinstance(self.left_end, int):
                if self.left_end > 0:
                    df = df[df['Base'] <= self.left_end] # subset
            # save data.frame
            df.to_json(self.plot_json)
            return df


    def fa_to_fq(self, x):
        pass


    def is_fasta(self, x):
        pass


    def is_fastq(self, x):
        pass


    def is_tool(self, name):
        return which(name) is not None


    def run_falco(self):
        """
        options:
        --outdir             Create all output files in the specified
        --nogroup            Disable grouping of bases for reads >50bp
        ## falco version before 1.2.1
        # --skip-html          Skip generating HTML file
        # --skip-short-summary Skip short summary
        ## falco version 1.2.1
        -skip-report         [Falco only] Do not create FastQC report HTML file.
        -skip-summary        [Falco only] Do not create FastQC summary file
        --quiet              Do not print more run info
        ## expect output
        out_dir/fastqc_data.txt
        """
        cmd = ' '.join([
            which('falco'),
            '--nogroup -skip-report -skip-summary --quiet',
            '--outdir {}'.format(self.out_dir),
            self.fq
        ])
        # out txt
        if os.path.exists(self.out_txt) and not self.overwrite:
            print('file exists: {}'.format(self.out_txt))
        else:
            try:
                run_shell_cmd(cmd)
                out_tmp = os.path.join(self.out_dir, 'fastqc_data.txt')
                if os.path.exists(out_tmp):
                    shutil.copy(out_tmp, self.out_txt)
                    os.remove(out_tmp)
            except:
                print('Failed to run Falco')


    def run_fastqc(self):
        """
        options:
        --outdir Create all output files in the specified output directory.
        --nogroup
        --extract
        --quiet
        ## expect output
        {filename}_fastqc/fastqc_data.txt
        """
        cmd = ' '.join([
            which('fastqc'),
            '--nogroup --extract --quiet',
            '--outdir {}'.format(self.out_dir),
            self.fq
        ])
        # out txt
        if os.path.exists(self.out_txt) and not self.overwrite:
            print('file exists: {}'.format(self.out_txt))
        else:
            try:
                run_shell_cmd(cmd)
                tmp_dir = os.path.join(self.out_dir, self.name+'_fastqc')
                out_tmp = os.path.join(tmp_dir, 'fastqc_data.txt')
                if os.path.exists(out_tmp):
                    shutil.copy(out_tmp, self.out_txt)
                    shutil.rmtree(tmp_dir)
                    # remove html and zip
                    if os.path.exists(tmp_dir+'.html'):
                        os.remove(tmp_dir+'.html')
                    if os.path.exists(tmp_dir+'.zip'):
                        os.remove(tmp_dir+'.zip')
            except:
                print('Failed to run FastQC')


    def pick_tool(self):
        if self.is_tool('falco'):
            out = self.run_falco
        elif self.is_tool('fastqc'):
            out = self.run_fastqc
        else:
            msg = '\n'.join(
                '='*80,
                'Required tool not found',
                'Try to install either of falco or fastqc',
                'see: https://github.com/smithlabcode/falco',
                'or: https://github.com/s-andrews/FastQC',
                'in brief:',
                '$ conda install -c bioconda falco'
                'or',
                '$ conda install -c bioconda fastqc',
                '='*80,
            )
            print(msg)
            out = None
        return out


    def plot(self):
        if os.path.exists(self.out_png) and not self.overwrite:
            print('file exists: {}'.format(self.out_png))
            return None
        try:
            # data
            df = self.get_df()
            # fig
            fig = px.bar(
                df, x='Base', y=['A', 'C', 'G', 'T'],
                color_discrete_sequence=['#11dd11', '#0c0cde', '#414141', '#e12626'],
                title='Per Base Content'
                )
            fig.update_layout(
                title="Per Base Content",
                xaxis_title="Position in sequence",
                yaxis_title="Percentage %",
                legend_title="Base",
                font=dict(
                    family="Courier New, monospace",
                    size=12,
                    color='#1B262C',
                    # color="RebeccaPurple"
                )
            )
            fig.write_image(self.out_png)
            print('Save to file: {}'.format(self.out_png))
        except:
            print('Failed to plot: {}'.format(self.out_png))


    def run(self):
        qc = self.pick_tool()
        if qc:
            qc()
            self.plot()
        # remove
        if self.rm_tmp:
            for i in [self.out_txt, self.out_json]: # , self.plot_json]:
                if os.path.exists(i):
                    os.remove(i)


class BaseContent(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        # self.fq = fq
        self.init_args()


    def init_args(self):
        args = {
            'fq': None,
            'out_dir': None,
            'left_end': 0,
            'parallel_jobs': 1,
            'overwrite': False,
            'rm_tmp': False,
        }
        self = update_obj(self, args, force=False)
        self.out_dir = fix_out_dir(self.out_dir)
        self.init_fq()


    def init_fq(self):
        if isinstance(self.fq, str):
            self.fq = [self.fq]
        # check if fq exists
        if isinstance(self.fq, list):
            self.fq = [i for i in self.fq if self.is_fastq(i)]


    def is_fastq(self, x):
        if isinstance(x, str):
            p = re.compile('.f(ast)?q(.gz)?$', flags=re.IGNORECASE)
            return isinstance(p.search(x), re.Match)
        elif isinstance(x, list):
            # return [self.is_fastq(i) for i in x]
            return list(map(self.is_fastq, x))
        else:
            return False


    def run_single_fx(self, i):
        args = self.__dict__.copy()
        args.pop('fq', []) # remove fq from args
        BaseContentR1(i, **args).run()


    def run(self):
        if len(self.fq) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fx, self.fq)
        else:
            [self.run_single_fx(i) for i in self.fq]


def base_content(**kwargs):
    # default args
    args = {
        'fq': None,
        'out_dir': None,
        'left_end': 0,
        'parallel_jobs': 1,
        'overwrite': False,
        'rm_tmp': False,
    }
    args.update(kwargs)
    fq = args.pop('fq', [])
    # BEGIN: sub-func #
    def run_single_fx(i):
        BaseContent(fq=i, **args).run()
    # END: sub-func #
    if len(fq) > 1 and args['parallel_jobs'] > 1:
        with Pool(processes=args['parallel_jobs']) as pool:
            pool.map(run_single_fx, fq)
    else:
        # print(kwargs)
        [BaseContent(fq=i, **args).run() for i in fq]


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. check base-content for single file',
        '$ python base_content.py -i fq1 -o out_dir',
        '2. check base-content for left 10-nt',
        '$ python base_content.py -i fq1 -o out_dir -l 10',
    ])
    parser = argparse.ArgumentParser(
        prog='base_content',
        description='Per base content for fastq file',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--fq', dest='fq', nargs='+', required=True,
        help='FASTQ files')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=False,
        help='The directory to save output files')
    parser.add_argument('-l', '--left-end', dest='left_end', type=int, default=0,
        help='Plot for the left N-base, default: [0], for whole sequence')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs', type=int,
        default=1, help='Run N-jobs in parallel, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite the exists files')
    parser.add_argument('-r', '--rm-tmp', dest='rm_tmp', action='store_true',
        help='remove temp files')
    return parser


def main():
    args = vars(get_args().parse_args())
    BaseContent(**args).run()


if __name__ == '__main__':
    main()

#