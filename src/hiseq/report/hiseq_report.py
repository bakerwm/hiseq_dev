#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
# import importlib
import argparse
import shutil
import importlib.resources as res
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd #, list_pkg_file
from hiseq.utils.file import file_exists, check_dir
from hiseq.utils.hiseq_utils import read_hiseq


class HiSeqRpt(object):
    """
    HiSeq reporter:
    Parameters
    ----------
    x:  str
        Path to the hiseq project_dir
    >>> HiSeqRpt(project_dir).run()
    >>> HiSeqRpt(project_dir).report()
    """
    def __init__(self, x, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.x = x
        self.overwrite = kwargs.get('overwrite', False)
        self.init_args()
        self.project_dir = x


    def init_args(self):
        if not file_exists(self.x):
            raise ValueError('project_dir not valid: {}'.format(self.x))
        # check hiseq_dir
        a = read_hiseq(self.x)
        if not a.is_hiseq:
            raise ValueError('project_dir not hiseq_dir: {}'.format(self.x))
        # check default files
        self.hiseq_type = 'trim_rp'
        self.report_dir = os.path.join(self.x, 'report')
        self.config_yaml = os.path.join(self.report_dir, 'config.yaml')
        self.report_html = os.path.join(self.report_dir, 'HiSeq_report.html')
        self.report_stdout = os.path.join(self.report_dir, 'report.stdout')
        self.report_stderr = os.path.join(self.report_dir, 'report.stderr')
        check_dir(self.report_dir)
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def report(self):
        hiseq_report_r = list_pkg_file('hiseq.data', 'hiseq_report.R')
        if not file_exists(hiseq_report_r):
            log.error('file not found: {}'.format(hiseq_report_r))
            return None
        # build command line
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            hiseq_report_r,
            self.project_dir,
            self.report_dir,
            '1>{}'.format(self.report_stdout),
            '2>{}'.format(self.report_stderr),
            ])
        # save command
        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # report_html
        if file_exists(self.report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                self.report_html))
        else:
            run_shell_cmd(cmd)
        # check again
        if not file_exists(self.report_html):
            log.error('report() failed, check log: {}'.format(self.report_stderr))


    def run(self):
        self.report()


def list_pkg_file(*args):
    """
    list the file in package
    >>> list_pkg_file('hiseq.data', 'illumina_index.csv')
    see importlib.resources.path()
    """
    try:
        with res.path(*args) as f:
            out = str(f)
    except:
        print('could not find file: {}'.format(args[-1]))
        out = None
    return out


def get_args():
    example = '\n'.join([
        'Examples:',
        '$ python hiseq_report.py -i project_dir',
    ])
    parser = argparse.ArgumentParser(
        prog='hiseq_report',
        description='Generate hiseq report',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--project-dir', dest='project_dir',
        required=True,
        help='Directory of hiseq project')
    parser.add_argument('-r', '--overwrite', action='store_true',
        help='Overwrite the exists file.')
    return parser


def main():
    args = get_args().parse_args()
    HiSeqRpt(x=args.project_dir, overwrite=args.overwrite).run()


if __name__ == '__main__':
    main()

#