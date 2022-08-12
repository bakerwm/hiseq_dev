#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download files from different sources:
1. OSS (aliyun, ossutil64) 
2. OSS (novogene, linuxnd, lnd)
2. url
"""

import os
# import sys
import re
# from this import d
# import yaml
import pathlib
import argparse
import subprocess
from urllib import request
# from datetime import datetime
from dateutil import tz
from hiseq.utils.utils import log, update_obj, get_date, Config
from hiseq.utils.file import file_abspath, check_dir


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
            'oss_txt': None,
            'out_dir': None,
            'region': 'shanghai',
            'ossutil': '/data/biosoft/ossutil_aliyuncs/ossutil64',
        }
        self = update_obj(self, args, force=False)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.url_list = self.parse_url(self.oss_txt)
        self.endpoint = 'http://oss-cn-{}.aliyuncs.com'.format(self.region.lower())
    

    def parse_url(self, x):
        url_list = []
        try:
            with open(x) as r:
                for l in r:
                    if l.startswith('http'):
                        url_list.append(l.strip())
        except:
            print('failed reading file: {}'.format(x))
        return url_list


    def get_url_name(self, x):
        """
        get the filename from url
        http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
        re.sub('^\w+\%|\?.*$', '', b)
        """
        if isinstance(x, str):
            return re.sub('^\w+\%|\?.*$', '', os.path.basename(x))


    def download_url(self, url, filename=None):
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
        if not isinstance(filename, str):
            filename = self.get_url_name(url)
        local_file = os.path.join(self.out_dir, filename)
        s = isinstance(url, str) and isinstance(local_file, str)
        if not s:
            return None
        # show msg
        msg = '\n'.join([
            '='*80,
            '{:>14} : {}'.format('Program', 'download_url'),
            '{:>14} : {}'.format('Date', get_date()),
            '{:>14} : {}'.format('url', url),
            '{:>14} : {}'.format('local file', local_file),
            '='*80,
        ])
        print(msg)
        # # config
        # local_dir = os.path.dirname(local_file)
        # if not os.path.exists(local_dir):
        #     os.makedirs(local_dir)
        check_dir(self.out_dir)
        if os.path.exists(local_file):
            print('local file exists')
        else:
            try:
                request.urlretrieve(url, local_file)
            except:
                print('failed downloading file: {}'.format(url))
        if not os.path.exists(local_file):
            print('failed downloading file: {}'.format(local_file))
        return local_file


    def parse_oss(self, x):
        """
        parse for AccessKeyId, AccessKeySecret, 预设OSS路径
        template:
        合同号	YF...
        开题单号	LPL...
        提取码	13a6
        AccessKeyId	...
        AccessKeySecret	...
        预设OSS路径	oss://...
        过期时间	2022/08/03 18:20:51
        使用教程	在线教程
        """
        # CN to EN
        k = {
            '合同号': 'contract_id',
            '开题单号': 'project_id',
            'AccessKeyId': 'AccessKeyId',
            'AccessKeySecret': 'AccessKeySecret',
            '预设OSS路径': 'OSS_path',
            '过期时间': 'expire',
        }
        # d = {'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
        d = {}
        try:
            with open(x) as r:
                for l in r:
                    s = l.strip().split(None, 1)
                    if len(s) == 2:
                        if s[0] in k:
                            kk = k.get(s[0])
                            d.update({kk:s[1]})
        except:
            print('failed reading file: {}'.format(x))
            d = None
        # # save as yaml
        # try:
        #     cfg = os.path.join(out_dir, 'config.yaml')
        #     if not os.path.exists(out_dir):
        #         os.makedirs(out_dir)
        #     with open(cfg, 'wt') as w:
        #         yaml.dump(d, w)
        # except:
        #     print('write to file, failed: {}'.format(cfg))
        return d


    # use ossutil command-line tool
    def list_oss(self, x):
        """
        List OSS files using ossutil64
        Parameters:
        -----------
        x : str
            text file, save the email content, see parse_oss()
        """
        print('!A-3', x)
        oss_cfg = self.parse_oss(x)
        cmd = ' '.join([
            self.ossutil,
            'ls',
            '-e {}'.format(self.endpoint),
            '-i {}'.format(oss_cfg.get('AccessKeyId')),
            '-k {}'.format(oss_cfg.get('AccessKeySecret')),
            oss_cfg.get('OSS_path'),
        ])
        # 3. run
        cmd_list = cmd.split()
        result = subprocess.run(cmd_list, stdout=subprocess.PIPE)
        # 4. check output
        if result.returncode == 0:
            s = result.stdout.decode('utf-8')
            # check fq.gz files #
            s1 = [i for i in s.split('\n') if i.endswith('.fq.gz')] # list rec
            s2 = [i.split()[-1] for i in s1] # file path
            s3 = list(map(os.path.basename, s2)) # filename
            s3a = [i for i in s3 if i.endswith('1.fq.gz')] # read1
            # s3b = [i for i in s3 if i.endswith('2.fq.gz')] # read2
            s3a.append('Number of fastq files (R1): {}'.format(len(s3a)))
            msg = '\n'.join(s3a)
            print(msg)
        else:
            ex = oss_cfg.get('expire', 'null')
            print('Could not access OSS: {}; expired at: {}'.format(x, ex))
        return result.returncode


    # use ossutil command-line tool
    def download_oss(self, x):
        """
        Download OSS files using ossutil64
        Parameters:
        -----------
        x : str
            text file, save the email content, see parse_oss()
        """
        s = self.list_oss(x) #
        if s > 0:
            print('error, could not list files: {}'.format(x))
            return None
        check_dir(self.out_dir)
        # 1. config
        oss_cfg = self.parse_oss(x)
        config_yaml = os.path.join(self.out_dir, 'config.yaml')
        Config().dump(oss_cfg, config_yaml)
        # 2. build command
        cmd = ' '.join([
            self.ossutil,
            'sync',
            '-e {}'.format(self.endpoint),
            '-i {}'.format(oss_cfg.get('AccessKeyId')),
            '-k {}'.format(oss_cfg.get('AccessKeySecret')),
            oss_cfg.get('OSS_path'),
            self.out_dir
        ])
        cmd_txt = os.path.join(self.out_dir, 'run.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # 3. run    
        msg = '\n'.join([
            '='*80,
            '{:>14} : {}'.format('program', 'download_oss'),
            '{:>14} : {}'.format('date', get_date()),
            '{:>14} : {}'.format('oss_path', oss_cfg.get('OSS_path')),
            '{:>14} : {}'.format('expire', oss_cfg.get('expire')),
            '='*80,
        ])
        print(msg)
        cmd_list = cmd.split()
        result = subprocess.run(cmd_list, stderr=subprocess.PIPE)
        if result.returncode > 0:
            ex = oss_cfg.get('expire', 'null')
            print('Could not access OSS: {}; expired at: {}'.format(x, ex))


    def run(self):
        url_list = self.parse_url(self.oss_txt)
        if len(url_list) > 0:
            [self.download_url(url) for url in url_list]
        else:
            self.download_oss(self.oss_txt)


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
    parser.add_argument('-i', '--oss-txt', dest='oss_txt',
        help='oss.txt file, could be a URL, or OSS config file')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
        help='directory to save the results')
    parser.add_argument('-r', '--region', dest='region', 
        default='shanghai',
        help='choose    region for public endpoint, default [shanghai]')
    parser.add_argument('--ossutil', 
        default='/data/biosoft/ossutil_aliyuncs/ossutil64',
        help='path to the ossutil command-line tool, default: [/data/biosoft/ossutil_aliyuncs/ossutil64]')
    return parser


def main():
    args = vars(get_args().parse_args())
    # download(**args)
    Download(**args).run()


if __name__ == '__main__':
    main()

# EOF

# # parse config from txt file
# def parse_oss(x, out_dir=None):
#     """
#     parse for AccessKeyId, AccessKeySecret, 预设OSS路径
#     template:
#     合同号	YF...
#     开题单号	LPL...
#     提取码	13a6
#     AccessKeyId	...
#     AccessKeySecret	...
#     预设OSS路径	oss://...
#     过期时间	2022/08/03 18:20:51
#     使用教程	在线教程
#     """
#     # CN to EN
#     k = {
#         '合同号': 'contract_id',
#         '开题单号': 'project_id',
#         'AccessKeyId': 'AccessKeyId',
#         'AccessKeySecret': 'AccessKeySecret',
#         '预设OSS路径': 'OSS_path',
#         '过期时间': 'expire',
#     }
#     d = {'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
#     try:
#         with open(x) as r:
#             for l in r:
#                 s = l.strip().split(None, 1)
#                 if len(s) == 2:
#                     if s[0] in k:
#                         kk = k.get(s[0])
#                         d.update({kk:s[1]})
#     except:
#         print('failed reading file: {}'.format(x))
#         d = None
#     # save as yaml
#     try:
#         cfg = os.path.join(out_dir, 'config.yaml')
#         if not os.path.exists(out_dir):
#             os.makedirs(out_dir)
#         with open(cfg, 'wt') as w:
#             yaml.dump(d, w)
#     except:
#         print('write to file, failed: {}'.format(cfg))
#     return d


# def get_date(timestamp=False):
#     now = datetime.now(tz.tzlocal())
#     if isinstance(timestamp, bool) and timestamp:
#         out = now.timestamp()
#     else:
#         now = now.astimezone(tz.tzlocal()) # to local
#         out = now.strftime('%Y-%m-%d %H:%M:%S') # YY-mm-dd H:M:S
#     return out


# # use ossutil commandline tool
# def list_oss(x, out_dir):
#     """
#     List OSS files using ossutil64
#     Parameters:
#     -----------
#     x : str
#         text file, save the email content, see parse_oss()
#     """
#     # 1. load config files
#     x = os.path.abspath(x)
#     out_dir = os.path.abspath(out_dir)
#     cfg = os.path.join(os.path.dirname(x), 'config.yaml')
#     d = parse_oss(x, out_dir)
#     # print(x, c, d)
#     # 2. prep command
#     ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil64' # !!! to-do
#     cmd = ' '.join([
#         ossutil,
#         'ls',
#         '-e {}'.format(d.get('endpoint')),
#         '-i {}'.format(d.get('AccessKeyId')),
#         '-k {}'.format(d.get('AccessKeySecret')),
#         d.get('OSS_path'),
#     ])
#     # 3. run
#     cmd_list = cmd.split()
#     result = subprocess.run(cmd_list, stdout=subprocess.PIPE)
#     if result.returncode == 0:
#         s = result.stdout.decode('utf-8')
#         # check fq.gz files #
#         s1 = [i for i in s.split('\n') if i.endswith('.fq.gz')] # list rec
#         s2 = [i.split()[-1] for i in s1] # file path
#         s3 = list(map(os.path.basename, s2)) # filename
#         s3a = [i for i in s3 if i.endswith('1.fq.gz')] # read1
#         # s3b = [i for i in s3 if i.endswith('2.fq.gz')] # read2
#         s3a.append('Number of fastq files (R1): {}'.format(len(s3a)))
#         msg = '\n'.join(s3a)
#         print(msg)
#     else:
#         ex = d.get('expire', 'null')
#         print('Could not access OSS: {}; expired at: {}'.format(x, ex))
#     return result.returncode


# # use ossutil command-line tool
# def download_oss(x, out_dir):
#     """
#     Download OSS files using ossutil64
#     Parameters:
#     -----------
#     x : str
#         text file, save the email content, see parse_oss()
#     """
#     # list files
#     s = list_oss(x, out_dir)
#     if s > 1:
#         print('error, could not list files: {}'.format(x))
#         return None
#         # download_oss(x, out_dir)
#     # 1. load config files
#     x = os.path.abspath(x)
#     out_dir = os.path.abspath(out_dir)
#     cfg = os.path.join(os.path.dirname(x), 'config.yaml')
#     d = parse_oss(x, out_dir)
#     # print(x, c, d)
#     # 2. prep command
#     ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil64' # !!! to-do
#     cmd = ' '.join([
#         ossutil,
#         'sync',
#         '-e {}'.format(d.get('endpoint')),
#         '-i {}'.format(d.get('AccessKeyId')),
#         '-k {}'.format(d.get('AccessKeySecret')),
#         d.get('OSS_path'),
#         out_dir
#     ])
#     cmd_txt = os.path.join(out_dir, 'run.sh')
#     with open(cmd_txt, 'wt') as w:
#         w.write(cmd+'\n')
#     # 3. run    
#     msg = '\n'.join([
#         '='*80,
#         '{:>14} : {}'.format('program', 'download_oss'),
#         '{:>14} : {}'.format('date', get_date()),
#         '{:>14} : {}'.format('oss_path', d.get('OSS_path')),
#         '{:>14} : {}'.format('expire', d.get('expire')),
#         '='*80,
#     ])
#     print(msg)
#     cmd_list = cmd.split()
#     result = subprocess.run(cmd_list, stderr=subprocess.PIPE)
#     if result.returncode > 0:
#         ex = d.get('expire', 'null')
#         print('Could not access OSS: {}; expired at: {}'.format(x, ex))


# def parse_url(x):
#     # parse oss.txt
#     url_list = []
#     try:
#         with open(x) as r:
#             for l in r:
#                 if l.startswith('http'):
#                     url_list.append(l.strip())
#                     # url = l.strip()
#                     # download_url(url, out_dir)
#     except:
#         print('failed reading file: {}'.format(x))
#     return url_list


# def get_url_name(x):
#     """
#     get the filename from url
#     http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
#     re.sub('^\w+\%|\?.*$', '', b)
#     """
#     if isinstance(x, str):
#         return re.sub('^\w+\%|\?.*$', '', os.path.basename(x))


# def download_url(url, out_dir, filename=None):
#     """
#     Download url and save to file    
#     from urllib import request
#     # Define the remote file to retrieve
#     remote_url = 'https://www.google.com/robots.txt'
#     # Define the local filename to save data
#     local_file = 'local_copy.txt'
#     # Download remote and save locally
#     request.urlretrieve(remote_url, local_file)
#     """    
#     if not isinstance(out_dir, str):
#         out_dir = str(pathlib.Path.cwd())
#     out_dir = os.path.abspath(os.path.expandvars(out_dir))
#     if not isinstance(filename, str):
#         filename = get_url_name(url)
#     local_file = os.path.join(out_dir, filename)
#     s = isinstance(url, str) and isinstance(local_file, str)
#     if not s:
#         return None
#     # show msg
#     msg = '\n'.join([
#         '='*80,
#         '{:>14} : {}'.format('program', 'download_url'),
#         '{:>14} : {}'.format('date', get_date()),
#         '{:>14} : {}'.format('url', url),
#         '{:>14} : {}'.format('local file', local_file),
#         '='*80,
#     ])
#     print(msg)
#     # config
#     local_dir = os.path.dirname(local_file)
#     if not os.path.exists(local_dir):
#         os.makedirs(local_dir)
#     if os.path.exists(local_file):
#         print('local file exists')
#     else:
#         try:
#             request.urlretrieve(url, local_file)
#         except:
#             print('failed downloading file: {}'.format(url))
#     if not os.path.exists(local_file):
#         print('failed downloading file: {}'.format(local_file))


# def download(oss_txt, out_dir):
#     url_list = parse_url(oss_txt)
#     if len(url_list) > 0:
#         [download_url(url, out_dir+'/'+get_url_name(url)) for url in url_list]
#     else:
#         download_oss(oss_txt, out_dir)

