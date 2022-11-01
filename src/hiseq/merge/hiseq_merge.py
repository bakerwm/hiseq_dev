#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Combine/merge multiple cnr report
cnr_r1, cnr_rn, ...
"""

import os
import pathlib
import argparse
# from hiseq.cnr.cnr_rp import CnrRp
from hiseq.utils.file import check_dir, file_abspath, list_dir, fix_out_dir
from hiseq.utils.utils import log, update_obj, Config
from hiseq.utils.hiseq_utils import is_hiseq_dir
from hiseq.cnr.cnr_utils import copy_hiseq_qc
from hiseq.atac.atac_files import get_atac_files
from hiseq.report.hiseq_report import HiSeqRpt


class HiSeqMerge(object):
    def __init__(self, **kwargs):
        c = HiSeqMergeConfig(**kwargs)
        self = update_obj(self, c.__dict__, force=True)
        
        
    def run(self):
        Config().dump(self.__dict__, self.config_yaml) # save config
        Config().dump({'in_dir': self.in_dir}, self.in_dir_json) # save dir list
        copy_hiseq_qc(self.project_dir) # copy files, config
        HiSeqRpt(self.project_dir).run() # generate report


class HiSeqMergeConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'in_dir': None, # required
            'out_dir': None,
            'merge_type': 'auto',
            'smp_name': 'hiseq_merge',
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_merge'
        self.out_dir = fix_out_dir(self.out_dir)
        self.in_dir = file_abspath(self.in_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.hiseq_list = self.load_input(self.in_dir, self.hiseq_type)
        self.init_files()


    def load_input(self, x, merge_type='auto'):
        """
        load hiseq dirs from input: [str, txt_file, json_file"hiseq_list"] ...
        return hiseq_dirs
        """
        out = []
        if isinstance(x, str):
            if os.path.isdir(x):
                out.append(x)
            elif os.path.isfile(x):
                if any([x.endswith(i) for i in ['.json', '.yaml', '.toml']]):
                    dx = Config().load(x)
                    out.extend(dx.get('hiseq_list', []))
                elif x.endswith('.txt'):
                    d = []
                    with open(x) as r:
                        for line in r:
                            line = line.strip()
                            if line.startswith('#') or len(line) == 0:
                                continue
                            d.append(line.split()[0]) # 
                    out.extend(d)
                else:
                    log.error(f'unknown input, expect [.txt, .json, .yaml, .toml], got {x}')
        elif isinstance(x, list):
            out = [j for i in x for j in self.load_input(i, merge_type)]
        else:
            log.error(f'unknown input, expect str, list, got {type(x).__name__}')
        # input should be hiseq_dir
        return [i for i in out if is_hiseq_dir(i, 'auto')]


    def init_files(self):
        # dirs
        self.project_dir = os.path.join(self.out_dir, self.smp_name)
        default_dirs = {
            'config_dir': 'config',
            'data_dir': 'data',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        ff = get_atac_files(self.out_dir, self.smp_name, None, None)
        fx = [
            'config_yaml', 'report_log', 'report_html', 'trim_summary_json',
            'align_summary_json', 'dup_summary_json', 'lendist_csv', 'lendist_txt',
            'lendist_pdf', 'frip_json'
        ]
        fd = {i:ff.get(i, None) for i in fx}
        fd.update({'in_dir_json': os.path.join(self.data_dir, 'in_dir.json')})
        self = update_obj(self, fd, force=True)
        # dirs
        dir_list = [
            self.config_dir, self.data_dir, self.qc_dir, self.report_dir,
        ]
        check_dir(dir_list)


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. Organize multiple hiseq dirs (cnr)',
        '$ python hiseq_merge.py -o ',
    ])
    parser = argparse.ArgumentParser(
        prog='hiseq_merge',
        description='hiseq_merge: merge multiple hiseq dirs',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--in-dir', dest='in_dir', nargs='+', required=True,
        help='hiseq project directories')
    parser.add_argument('-o', '--out-dir', dest='out_dir', default=None,
        help='Directory saving results, default: [./hiseq_merge]')
    parser.add_argument('-n', '--smp-name', dest='smp_name', default='hiseq_merge',
        help='Name of the hiseq_merge files, default: [hiseq_merge]')
    parser.add_argument('-f', '--force', dest='overwrite', action='store_true',
        help='ignore exists files')
    return parser


def main():
    args = vars(get_args().parse_args())
    HiSeqMerge(**args).run()


if __name__ == '__main__':
    main()

#
