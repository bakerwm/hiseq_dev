#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download files from different sources:
1. OSS (aliyun, ossutil64) 
2. OSS (novogene, linuxnd, lnd)
2. url
"""

import os
import sys
import re
import yaml
import argparse
import subprocess
from urllib import request
from datetime import datetime
from dateutil import tz
from hiseq.utils.utils import log


# parse config from txt file
def parse_oss(x, out_dir=None):
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
    d = {'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
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
    # save as yaml
    try:
        cfg = os.path.join(out_dir, 'config.yaml')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(cfg, 'wt') as w:
            yaml.dump(d, w)
    except:
        print('write to file, failed: {}'.format(cfg))
    return d


# use ossutil commandline tool
def download_oss(x, out_dir):
    """
    Download OSS files using ossutil64
    Parameters:
    -----------
    x : str
        text file, save the email content, see parse_oss()
    """
    # 1. load config files
    x = os.path.abspath(x)
    out_dir = os.path.abspath(out_dir)
    cfg = os.path.join(os.path.dirname(x), 'config.yaml')
    d = parse_oss(x, out_dir)
    # print(x, c, d)
    # 2. prep command
    ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil64' # !!! to-do
    cmd = ' '.join([
        ossutil,
        'sync',
        '-e {}'.format(d.get('endpoint')),
        '-i {}'.format(d.get('AccessKeyId')),
        '-k {}'.format(d.get('AccessKeySecret')),
        d.get('OSS_path'),
        out_dir
    ])
    cmd_txt = os.path.join(out_dir, 'run.sh')
    with open(cmd_txt, 'wt') as w:
        w.write(cmd+'\n')
    # 3. run    
    msg = '\n'.join([
        '='*80,
        '{:>14} : {}'.format('program', 'download_oss'),
        '{:>14} : {}'.format('date', get_date()),
        '{:>14} : {}'.format('oss_path', d.get('OSS_path')),
        '{:>14} : {}'.format('expire', d.get('expire')),
        '='*80,
    ])
    print(msg)
    cmd_list = cmd.split()
    result = subprocess.run(cmd_list, stderr=subprocess.PIPE)
    if result.returncode > 0:
        ex = d.get('expire', 'null')
        print('Could not access OSS: {}; expired at: {}'.format(x, ex))


# use ossutil commandline tool
def list_oss(x, out_dir):
    """
    List OSS files using ossutil64
    Parameters:
    -----------
    x : str
        text file, save the email content, see parse_oss()
    """
    # 1. load config files
    x = os.path.abspath(x)
    out_dir = os.path.abspath(out_dir)
    cfg = os.path.join(os.path.dirname(x), 'config.yaml')
    d = parse_oss(x, out_dir)
    # print(x, c, d)
    # 2. prep command
    ossutil = '/data/biosoft/ossutil_aliyuncs/ossutil64' # !!! to-do
    cmd = ' '.join([
        ossutil,
        'ls',
        '-e {}'.format(d.get('endpoint')),
        '-i {}'.format(d.get('AccessKeyId')),
        '-k {}'.format(d.get('AccessKeySecret')),
        d.get('OSS_path'),
    ])
    # 3. run
    cmd_list = cmd.split()
    result = subprocess.run(cmd_list, stdout=subprocess.PIPE)
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
        ex = d.get('expire', 'null')
        print('Could not access OSS: {}; expired at: {}'.format(x, ex))
    return result.returncode


def download_url(url, local_file):
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
    s = isinstance(url, str) and isinstance(local_file, str)
    if not s:
        return None
    # show msg
    msg = '\n'.join([
        '='*80,
        '{:>14} : {}'.format('program', 'download_url'),
        '{:>14} : {}'.format('date', get_date()),
        '{:>14} : {}'.format('url', url),
        '{:>14} : {}'.format('local file', local_file),
        '='*80,
    ])
    print(msg)
    # config
    local_dir = os.path.dirname(local_file)
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)
    if os.path.exists(local_file):
        print('local file exists')
    else:
        try:
            request.urlretrieve(url, local_file)
        except:
            print('failed downloading file: {}'.format(url))
    if not os.path.exists(local_file):
        print('failed downloading file: {}'.format(local_file))


def get_date(timestamp=False):
    now = datetime.now(tz.tzlocal())
    if isinstance(timestamp, bool) and timestamp:
        out = now.timestamp()
    else:
        now = now.astimezone(tz.tzlocal()) # to local
        out = now.strftime('%Y-%m-%d %H:%M:%S') # YY-mm-dd H:M:S
    return out


def get_url_name(x):
    """
    get the filename from url
    http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
    re.sub('^\w+\%|\?.*$', '', b)
    """
    if isinstance(x, str):
        return re.sub('^\w+\%|\?.*$', '', os.path.basename(x))


def download(oss_txt, out_dir):
    # parse oss.txt
    url_list = []
    with open(oss_txt) as r:
        for l in r:
            if l.startswith('http'):
                url_list.append(l.strip())
                # url = l.strip()
                # download_url(url, out_dir)
    # choose downloader
    if len(url_list) > 0:
        [download_url(url, out_dir+'/'+get_url_name(url)) for url in url_list]
    else:
        s = list_oss(sys.argv[1], sys.argv[2])
        if s == 0:
            download_oss(sys.argv[1], sys.argv[2])



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
    return parser


def main():
    # if len(sys.argv) < 3:
    #     print('Usage: python download.py <oss.txt> <out_dir>')
    #     sys.exit(1)
    # oss_txt, out_dir = sys.argv[1:3]
    args = vars(get_args().parse_args())
    download(**args)


if __name__ == '__main__':
    main()

# EOF