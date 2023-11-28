#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download files from different sources:
1. OSS (aliyun, ossutil64) 
2. OSS (novogene, linuxnd, lnd)
2. url
"""

import os
import re
import argparse
import subprocess
from pathlib import Path
from urllib import request
# from dateutil import tz
from hiseq.utils.utils import log, update_obj, get_date


class OSSinfo(object):
    """
    Parse OSS info from text file, that could be in the following format:
    format-1:
        AccessKeyId: ...
        Access_key_secret: ...
        预设OSS路径: oss://...
        区域: 华北2(北京)
    format-2:
        AccessKeyId     ...
        AccessKeySecret ...
        预设OSS路径     oss://...
    format-3:
        AccessKeyId ...
        AccessKeySecret ...
        预设TOS路径: tos://sky...
        Endpoint: https://tos-cn-shanghai.volces.com
        华东2(上海)
    format-4:
        https://...
    format-5:
        oss://...
        tos://...
    # return 
    1. list of oss_info
    2. list of http urls
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args = {
            'oss_txt': None,
            'region': None,
            'endpoint': None,
            'ossutil': '/data/biosoft/ossutil/ossutil',
            'accesskeyid': None,
            'accesskeysecret': None,
        }
        self = update_obj(self, args, force=False)
        self.get_ossutil() #
        self.info = [self.parse_oss_info(i) for i in self.load_oss_txt()]


    def get_ossutil(self):
        if isinstance(self.ossutil, str):
            self.tosutil = str(Path(self.ossutil).parent / 'tosutil')
        else:
            self.tosutil = 'tosutil' # default
        if not Path(self.ossutil).exists():
            self.ossutil = 'ossutil' # from PATH
        if not Path(self.tosutil).exists():
            self.tosutil = 'tosutil' # from PATH


    def load_oss_txt(self, x=None):
        """
        Load oss_txt file as list
        split by empty lines
        """
        if x is None:
            x = self.oss_txt
        blocks = [] # init
        try:
            with open(x) as r:
                line = r.read() # read whole file at once
                # split into blocks by empty line
                line = re.sub('\\n\\s*\\n', '\\n\\n', line.strip())
                blocks = re.split('\\n{2,100}', line)
                blocks = [i for i in blocks if len(i) > 20]
        except FileNotFoundError as e:
            log.warning(f'File not found: {e}, {x}')
        except Exception as e:
            log.warning(f'Could not read file: {e}, {x}')
        return blocks


    def parse_oss_info(self, x):
        """
        Parse OSS key_id, key_secret, ... from string
        return dict
        """
        # parse text
        if isinstance(x, str):
            # check version
            config = {
                'http': None,
                'accesskeyid': None,
                'accesskeysecret': None,
                'oss_path': None,
                'endpoint': None,
                'region': None,
            }
            for line in re.split('\\n', x.strip()):
                line = line.strip()
                if line.startswith('#') or len(line) < 20:
                    continue
                # single value
                if line.startswith('http://'):
                    config.update({'http':line})
                if line.startswith('tos://'):
                    config.update({'oss_path': line})
                # key:value pairs
                tabs = re.split('[\\t\\s=：]', line.strip(), maxsplit=1) # Chinese characters
                if len(tabs) < 2:
                    continue
                k,v = tabs[:2]
                k = k.lower().replace('_', '') # format
                v = v.strip() 
                # skip default empty values
                if '*' in v:
                    v = ''
                # guess attributes
                if v.startswith('oss://'):
                    config.update({
                        'oss_path': v, 
                        'oss_type': 'oss', 
                        'oss_config': str(Path('~/.ossutilconfig').expanduser()),
                        'ossutil': self.ossutil,
                    })
                elif v.startswith('tos://'):
                    config.update({
                        'oss_path': v, 
                        'oss_type': 'tos', 
                        'oss_config': str(Path('~/.tosutilconfig').expanduser()),
                        'ossutil': self.tosutil,
                    })
                elif 'keyid' in k or 'ak' == k:
                    config.update({'accesskeyid': v})
                elif 'secret' in k or 'sk' == k:
                    config.update({'accesskeysecret': v})
                elif 'endpoint' in k:
                    config.update({'endpoint': v})
                elif 'region' in k:
                    config.update({'region': v})
                else:
                    pass
            # load config from $HOME for tos/oss
            # update endpoint, region
            defaults = {
                'psndata': {'region': 'shanghai', 'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'},
                'seekgene': {'region': 'beijing', 'endpoint': 'http://oss-cn-beijing.aliyuncs.com'},
                'skyseq': {'region': 'cn-shanghai', 'endpoint': 'https://tos-cn-shanghai.volces.com'},
            }
            # check oss source
            for k,v in defaults.items():
                oss_path = config.get('oss_path', None)
                if isinstance(oss_path, str):
                    if k in oss_path:
                        if config.get('region', None) is None:
                            config.update({'region': v.get('region', None)})
                        if config.get('endpoint', None) is None:
                            config.update({'endpoint': v.get('endpoint', None)})
            return config


class Download(object):
    """
    Download files from url or OSS (Aliyun)    
    see https://help.aliyun.com/document_detail/31837.html
    for more details about public endpoints for Aliyun
    shanghai http://oss-cn-shanghai.aliyuncs.com
    beijing http://oss-cn-beijing.aliyuncs.com
    ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args = {
            'out_dir': None,
            'http': None,
            'accesskeyid': None,
            'accesskeysecret': None,
            'oss_path': None,
            'endpoint': None,
            'region': None,
            'ossutil': None,
            'overwrite': False,
        }
        self = update_obj(self, args, force=False)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd() / 'from_illumina')
        self.out_dir = str(Path(self.out_dir).expanduser().absolute())
        # self.is_http = isinstance(self.http, str)
        if isinstance(self.http, str) or self.is_valid_oss():
            if not Path(self.out_dir).exists():
                try:
                    Path(self.out_dir).mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(f'Could not create dir, {e}, {self.out_dir}')


    def is_valid_oss(self):
        opts = {
            'accesskeyid': self.accesskeyid,
            'accesskeysecret': self.accesskeysecret,
            'oss_path': self.oss_path,
            'endpoint': self.endpoint,
            'region': self.region,
            'ossutil': self.ossutil,
        }
        out = all([isinstance(i, str) for i in list(opts.values())])
        if not out:
            msg = '\n'.join([f'{k:<15} = {v}' for k,v in opts.items()])
            print(msg)
        return out


    def get_url_name(self, x):
        """
        get the filename from url
        http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
        re.sub('^\w+\%|\?.*$', '', b)
        """
        if isinstance(x, str):
            xname = re.sub('?Expires=.*', '', Path(x).name) # fix
            return re.sub('[^A-Za-z0-9\\_\\-\\.]', '', xname)
            # return re.sub('^\w+\%|\?.*$', '', os.path.basename(x))


    def download_http(self):
        """
        Download url and save to file    
        from urllib import request
        # Define the remote file to retrieve
        remote_url = 'https://www.google.com/robots.txt'
        # Define the local filename to save data
        local_file = 'local_copy.txt'
        # Download remote and save locally
        request.urlretrieve(remote_url, local_file)
        """
        if not isinstance(self.http, str):
            print(f'skipped, http is not str, got {type(self.http)}')
            return None
        filename = self.get_url_name(self.http)
        dest_file = str(Path(self.out_dir) / filename)
        # show msg    
        msg = '\n'.join([
            '='*80,
            f'{"Program":>14} : {"download_http"}',
            f'{"Date":>14} : {get_date()}',
            f'{"url":>14}: {self.download_http}',
            f'{"out_dir":>14}: {self.out_dir}',
            '='*80,
        ])
        print(msg)
        # check_dir(self.out_dir)
        if Path(dest_file).exists() and not self.overwrite:
            print(f'file exists: {dest_file}')
        else:
            try:
                request.urlretrieve(url, dest_file)
            except Exception as e:
                print(f'failed downloading url, {e}, {self.http}')
        if not Path(dest_file).exists():
            print(f'file not found: {dest_file}')
        return dest_file


    def get_oss_cmd(self, subcommand='ls', options=''):
        if self.is_valid_oss():
            # check ossutil command
            if not Path(self.ossutil).exists():
                print(f'ossutil command not found: {self.ossutil}')
                print(f'find tosutil at: "https://www.volcengine.com/docs/6349/148777"')
                print(f'find ossutil at: "https://help.aliyun.com/zh/oss/developer-reference/oss-tools"')
                return None
            # tos require 'region'
            if self.oss_path.startswith('tos://'):
                opt_region = f'-re {self.region}'
            else:
                opt_region = ''
            # dest
            if subcommand in ['cp', 'sync']:
                opt_dest = self.out_dir
            else:
                opt_dest = ''
            # command
            cmd = ' '.join([
                self.ossutil,
                subcommand,
                options,
                opt_region,
                f'-e {self.endpoint}',
                f'-i {self.accesskeyid}',
                f'-k {self.accesskeysecret}',
                self.oss_path,
                opt_dest,
            ])
            return cmd


    # use ossutil command-line tool
    def list_oss(self):
        cmd = self.get_oss_cmd('ls')
        if not isinstance(cmd, str):
            print(f'invalid oss_info')
            return 1 # exit code
        # 3. run
        result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
        # 4. check output
        if result.returncode > 0:
            print(result.stdout.decode())
        else:
            s = result.stdout.decode('utf-8')
            fq_list = [Path(i).name for i in re.split('\\n', s) if i.endswith('gz')]
            msg_list = [f'{i+1:>02}. {k}' for i,k in enumerate(fq_list)]
            msg_list.append(f'[{len(fq_list)}] fastq files found')
            msg = '\n'.join(msg_list)
            print(msg)
        return result.returncode # 0=yes, 1=no


    # use ossutil command-line tool
    def download_oss(self, subcommand='cp'):
        """
        Download OSS files using ossutil/tosutil
        subcommand: cp, sync
        """
        if self.list_oss() > 0: # return code
            print(f'failed to list oss info: {self.oss_path}')
            return None
        # build command
        cmd = self.get_oss_cmd('cp', options='-r')
        cmd_txt = str(Path(self.out_dir) / 'run.sh')
        try:
            with open(cmd_txt, 'wt') as w:
                w.write(cmd+'\n')
        except Exception as e:
            print(f'failed, could not write to file: {cmd_txt}')
            return None
        # show log
        msg = '\n'.join([
            '='*80,
            f'{"program":>14} : {"download_oss"}',
            f'{"date":>14} : {get_date()}',
            f'{"oss_path":>14}: {self.oss_path}',
            f'{"out_dir":>14}: {self.out_dir}',
            '='*80,
        ])
        print(msg)
        # run
        try:
            result = subprocess.run(cmd.split(), stderr=subprocess.PIPE)
            if result.returncode > 0:
                print(f'failed, ossutil {subcommand}, {result.stderr}')
        except Exception as e:
            print(f'ossutil failed, {e}')


    def run(self):
        if isinstance(self.http, str):
            self.download_http()
        else:
            self.download_oss()


def download(oss_txt, out_dir, **kwargs):
    kwargs.update({'oss_txt': oss_txt, 'out_dir': out_dir})
    for oss_info in OSSinfo(**kwargs).info:
        oss_info.update({
            'out_dir': out_dir, 
            'overwrite': kwargs.get('overwrite', False),
        })
        # print('!A-1', oss_info)
        dn = Download(**oss_info)
        dn.run() 


def get_args():
    example = '\n'.join([
        'Download fq data from companies, OSS, ...',
        '1. download oss data',
        '$ python download.py -i oss.txt -o out_dir',
    ])
    parser = argparse.ArgumentParser(
        prog='download',
        description='download oss data',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--oss-txt', dest='oss_txt', required=True,
        help='oss.txt file, could be a URL, or OSS config file')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    parser.add_argument('-r', '--region', dest='region', default=None,
        help='choose region, eg: beijing, shanghai, ..., default [None]')
    parser.add_argument('-e', '--endpoint', dest="endpoint", default=None,
        help='specify the endpoint, eg: http://oss-cn-shanghai.aliyuncs.com, default [None]')
    parser.add_argument('-t', '--ossutil', dest='ossutil',
        default='/data/biosoft/ossutil/ossutil64',
        help='ossutil command-line tool, default: [/data/biosoft/ossutil/ossutil64]')
    parser.add_argument('-i', '--AccessKeyId', dest='accesskeyid', default=None,
        help='The AccessKeyId, default: [None]')
    parser.add_argument('-k', '--AccessKeySecret', dest='accesskeysecret', default=None,
        help='The AccessKeySecret, default: [None]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files')
    return parser


def main():
    args = vars(get_args().parse_args())
    oss_txt = args.pop('oss_txt', None)
    out_dir = args.pop('out_dir', None)
    download(oss_txt, out_dir, **args)


if __name__ == '__main__':
    main()

# EOF


# def parse_oss_info(self, x):
#     """
#     Parse OSS key_id, key_secret, ... from string
#     return dict
#     """
#     # check version
#     config = {
#         'http': None,
#         'accesskeyid': None,
#         'accesskeysecret': None,
#         'oss_path': None,
#         'endpoint': None,
#         'region': None,
#     }
#     # parse text
#     try:
#         with open(x) as r:
#             for line in r:
#                 if line.startswith('#') or len(line) < 4:
#                     continue
#                 if 'http://geneplus' in line:
#                     config.update({'http':line.strip()})
#                 tabs = re.split('[\\t\\s=：]', line.strip(), maxsplit=1) # chinese characters
#                 if len(tabs) < 2:
#                     continue
#                 k,v = tabs[:2]
#                 k = k.lower().replace('_', '')
#                 # skip default empty values
#                 if '*' in v:
#                     v = ''
#                 # guess attributes
#                 if 'oss://' in line:
#                     oss_config = str(Path('~/.ossutilconfig').expanduser())
#                     config.update({'oss_path': v, 'oss_type': 'oss', 'oss_config': oss_config})
#                 elif 'tos://' in line:
#                     oss_config = str(Path('~/.tosutilconfig').expanduser())
#                     config.update({'oss_path': v, 'oss_type': 'tos', 'oss_config': oss_config})
#                 elif 'keyid' in k or 'ak' == k:
#                     config.update({'accesskeyid': v})
#                 elif 'secret' in k or 'sk' == k:
#                     config.update({'accesskeysecret': v})
#                 elif 'endpoint' in k:
#                     config.update({'endpoint': v})
#                 elif 'region' in k:
#                     config.update({'region': v})
#                 else:
#                     pass
#     except FileNotFoundError as e:
#         log.warning(f'File not found: {e}, {x}')
#     except Exception as e:
#         log.warning(f'Could not read file: {e}, {x}')
#     # load config from $HOME for tos/oss
#     # update endpoint, region
#     defaults = {
#         'psndata': {'region': 'shanghai', 'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'},
#         'seekgene': {'region': 'beijing', 'endpoint': 'http://oss-cn-beijing.aliyuncs.com'},
#         'skyseq': {'region': 'cn-shanghai', 'endpoint': 'https://tos-cn-shanghai.volces.com'},
#     }
#     # check oss source
#     for k,v in defaults.items():
#         oss_path = config.get('oss_path', None)
#         if isinstance(oss_path, str):
#             if k in oss_path:
#                 if config.get('region', None) is None:
#                     config.update({'region': v.get('region', None)})
#                 if config.get('endpoint', None) is None:
#                     config.update({'endpoint': v.get('endpoint', None)})
#     # pass
#     return config


