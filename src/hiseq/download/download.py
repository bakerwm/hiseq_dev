#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download files from different sources:
1. OSS (aliyun, ossutil) 
2. OSS (novogene, linuxnd, lnd)
2. url
"""

import os
import re
import argparse
import subprocess
from pathlib import Path
from urllib import request
from dateutil import tz
from datetime import datetime

# from hiseq.utils.utils import update_obj, get_date

from hiseq.download.ossutil import TOS, OSS


def update_obj(obj, d, force=True, remove=False):
    """
    Update the object, by dict
    Args:
        obj (object): an object
        d (dict): a dict to update values in object
        force (bool): update exists attributes to object, default: True
        remove (bool): remove exists attributes from object, default: False
    Returns:
        (object) an updated object
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj


def get_date(timestamp=False):
    """
    Return the current date in UTC.timestamp or local-formated-string
    calculation in UTC
    Args:
        timestamp (bool): return timestamp
    Returns:
        current time

    Example:
    >>> from datetime import datetime
    >>> from dateutil import tz
    >>> get_date()
    '2021-05-18 17:08:53'

    >>> get_date(True)
    1621328957.280303

    # convert timestamp to local time
    >>> ts = get_date(True)
    >>> datetime.fromtimestamp(ts, tz.tzlocal())
    """
    now = datetime.now(tz.tzlocal())
    if isinstance(timestamp, bool) and timestamp:
        out = now.timestamp()
    else:
        now = now.astimezone(tz.tzlocal())  # to local
        out = now.strftime("%Y-%m-%d %H:%M:%S")  # YY-mm-dd H:M:S
    return out


class OSSinfo(object):
    """
    Parse OSS info from text file, that could be in the following format:

    Args:
        oss_txt (str): A text file record the oss info
        ossutil (str): Path to the ossutil command
        endpoint (str): endpoint for oss
        region (str): region for oss
        accesskeyid (str): AccessKeyId for oss
        accesskeysecret (str): AccessKeyId for oss

    Returns:
        (list) a list of oss_info in dict

    Description:
    format-1:
        AccessKeyId: ...
        Access_key_secret: ...
        预设OSS路径: oss://...
        区域: 华北2(北京)
    format-2:
        AccessKeyId     ...
        AccessKeySecret ...
        预设OSS路径     oss://...
    format-4:
        https://...
    format-5:
        oss://...
        tos://...
    # pre-defined config
    'qgq': {'region': 'shanghai', 'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
    'psndata': {'region': 'shanghai', 'endpoint': 'http://oss-cn-shanghai.aliyuncs.com'}
    'seekgene': {'region': 'beijing', 'endpoint': 'http://oss-cn-beijing.aliyuncs.com'}
    'skyseq': {'region': 'cn-shanghai', 'endpoint': 'https://tos-cn-shanghai.volces.com'}
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args = {
            "oss_txt": None,
            "region": None,
            "endpoint": None,
            "ossutil": None,
            "accesskeyid": None,
            "accesskeysecret": None,
        }
        self = update_obj(self, args, force=False)
        self.info = [self.parse_oss_info(i) for i in self.load_oss_txt()]

    def pre_config(self, oss_path=None):
        """
        pre-defined config for selected companies
        oss_path: 'oss://seekgene-release/2307015_1010/release/'
        """
        cfg = {
            "qgq": {
                "region": "shanghai",
                "endpoint": "http://oss-cn-shanghai.aliyuncs.com",
            },
            "psndata": {
                "region": "shanghai",
                "endpoint": "http://oss-cn-shanghai.aliyuncs.com",
            },
            "seekgene": {
                "region": "beijing",
                "endpoint": "http://oss-cn-beijing.aliyuncs.com",
            },
            "skyseq": {
                "region": "cn-shanghai",
                "endpoint": "https://tos-cn-shanghai.volces.com",
            },
        }
        if isinstance(oss_path, str):
            parts = Path(oss_path).parts
            cname = parts[1] if len(parts) > 1 else "psndata"  #
        else:
            cname = "psndata"
        cinfo = [v for k, v in cfg.items() if cname.startswith(k)]
        return cinfo[0] if len(cinfo) > 0 else {}

    def cmd_config(self, oss_path=None):
        if isinstance(self.ossutil, str):
            cmd = Path(self.ossutil).name
            cmd = re.sub("[^A-Za-z\\_\\-]", "", cmd)
            # auto switch between ossutil and tosutil
            if isinstance(oss_path, str):
                cmd = f"{oss_path[:3]}util"  # ossutil/tosutil
                # update
                self.ossutil = str(Path(self.ossutil).parent / cmd)
            cfg = str(Path("~").expanduser() / f".{cmd}config")
            return cfg

    def filter_oss(self, x):
        """
        Filter oss_info
        remove comment lines
        """
        if isinstance(x, str):
            lines = [
                i.strip()
                for i in re.split("\\n", x)
                if not i.startswith("#") and len(i) > 20
            ]
            return "\\n".join(lines)

    def load_oss_txt(self, x=None):
        """
        Load oss_txt file and split by blank line

        Args:
            x (str): Path to oss_txt file, default from self.oss_txt

        Returns:
            list: oss records
        """
        if x is None:
            x = self.oss_txt
        blocks = []  # init
        try:
            with open(x) as r:
                line = r.read()  # read whole file at once
                # split into blocks by empty line
                line = re.sub("\\n\\s*\\n", "\\n\\n", line.strip())
                blocks = re.split("\\n{2,100}", line)
                blocks = [i for i in blocks if len(self.filter_oss(i)) > 20]
        except FileNotFoundError as e:
            print(f"File not found: {e}, {x}")
        except Exception as e:
            print(f"Could not read file: {e}, {x}")
        return blocks

    def load_config(self, x=None):
        """
        Load ~/.ossutilconfig or ~/.tosutilconfig
        """
        cfg = {}  #
        try:
            with open(x) as r:
                for line in r:
                    tabs = re.split("=", line.strip(), maxsplit=1)
                    if len(tabs) < 2:
                        continue
                    if len(tabs[1]) < 1 or tabs[1].startswith("*"):
                        continue
                    # parse key
                    k = re.sub("[^A-Za-z]", "", tabs[0].lower())  #
                    if k == "ak":
                        k = "accesskeyid"
                    elif k == "sk":
                        k = "accesskeysecret"
                    cfg.update({k: tabs[1]})
        except Exception as e:
            print(f"Failed reading config: {e}, {x}")
        return cfg

    def parse_oss_info(self, x):
        """
        Parse OSS key_id, key_secret, ... from string

        Args:
            x (str): the oss info

        Returns
            dict: oss info as dict
        """
        # parse text
        if not isinstance(x, str):
            print(f"failed, parse_oss_info(x), expect str, got {type(x)}")
            return None
        # info
        config = {
            "http": None,
            "accesskeyid": None,
            "accesskeysecret": None,
            "oss_path": None,
            "endpoint": None,
            "region": None,
            # 'ossutil': self.ossutil,
        }
        # parse lines
        lines = [
            i.strip()
            for i in re.split("\\n", x)
            if not i.startswith("#") and len(i) > 20
        ]
        if len(lines) < 3:
            # http mode # caution multiple lines, only last line saved
            for line in lines:
                tabs = line.strip().split(maxsplit=1)
                v = tabs[1] if len(tabs) > 1 else line.strip()
                if v.startswith("http://"):
                    config.update({"http": line})
                elif v[2:].startswith("s://"):
                    config.update({"oss_path": line})
        else:
            # oss mode, in-case Chinese characters
            tabs = [re.split("[\\t\\s=:：]", i, maxsplit=1) for i in lines]
            tabs = [i for i in tabs if len(i) == 2]
            for k, v in tabs:
                # k = k.lower().replace('_', '') # format
                k = re.sub("[^A-Za-z]", "", k.lower())
                v = v.strip()
                if len(v) < 2:
                    continue
                if v[2:].startswith("s://"):
                    k = "oss_path"
                elif "keyid" in k or "ak" == k:
                    k = "accesskeyid"
                elif "keysecret" in k or "sk" == k:
                    k = "accesskeysecret"
                elif "endpoint" in k:
                    k = "endpoint"
                elif "region" in k:
                    k = "region"
                else:
                    pass
                config.update({k: v})
        # pre-defined region, endpoint
        pre_cfg = self.pre_config(config.get("oss_path", None))
        # update ossutil, ossutilconfig
        c_cfg = self.cmd_config(config.get("oss_path", None))
        cfg = {}  # config skipped !!! to-do
        # cfg = self.load_config(c_cfg) if Path(c_cfg).exists() else {}
        pre_cfg.update(cfg)  # all
        # update config
        for k, v in pre_cfg.items():
            if k in config and config.get(k, None) is None:
                config.update({k: v})
        # cmd updated, see self.cmd_config()
        config.update({"ossutil": self.ossutil})
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
            "out_dir": None,
            "http": None,
            "accesskeyid": None,
            "accesskeysecret": None,
            "oss_path": None,
            "endpoint": None,
            "region": None,
            "ossutil": None,
            "overwrite": False,
            "dry_run": False,
        }
        self = update_obj(self, args, force=False)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd() / "from_illumina")
        self.out_dir = str(Path(self.out_dir).expanduser().absolute())
        # self.is_http = isinstance(self.http, str)
        if isinstance(self.http, str) or self.is_valid_oss():
            if not Path(self.out_dir).exists():
                try:
                    Path(self.out_dir).mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(f"Could not create dir, {e}, {self.out_dir}")

    def is_valid_oss(self, verbose=True):
        opts = {
            "accesskeyid": self.accesskeyid,
            "accesskeysecret": self.accesskeysecret,
            "oss_path": self.oss_path,
            "endpoint": self.endpoint,
            "region": self.region,
            "ossutil": self.ossutil,
        }
        # update cmd
        out = all([isinstance(i, str) for i in list(opts.values())])
        if isinstance(self.ossutil, str):
            if not Path(self.ossutil).exists():
                out = False
                print(f"ossutil command not exists: {self.ossutil}")
        else:
            out = False
            print(f"ossutil missing, see Download(ossutil=)")
        if not out:
            msg = "\n".join([f"{k:<15} = {v}" for k, v in opts.items()])
            if verbose:
                print(msg)
        return out

    def get_url_name(self, x):
        """
        get the filename from url
        http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
        re.sub('^\w+\%|\?.*$', '', b)
        """
        if isinstance(x, str):
            xname = re.sub("?Expires=.*", "", Path(x).name)  # fix
            return re.sub("[^A-Za-z0-9\\_\\-\\.]", "", xname)
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
            print(f"skipped, http is not str, got {type(self.http)}")
            return None
        filename = self.get_url_name(self.http)
        dest_file = str(Path(self.out_dir) / filename)
        # show msg
        msg = "\n".join(
            [
                "=" * 80,
                f'{"Program":>14} : {"download_http"}',
                f'{"Date":>14} : {get_date()}',
                f'{"url":>14}: {self.download_http}',
                f'{"out_dir":>14}: {self.out_dir}',
                "=" * 80,
            ]
        )
        print(msg)
        # check_dir(self.out_dir)
        if Path(dest_file).exists() and not self.overwrite:
            print(f"file exists: {dest_file}")
        else:
            try:
                request.urlretrieve(self.http, dest_file)
            except Exception as e:
                print(f"failed downloading url, {e}, {self.http}")
        if not Path(dest_file).exists():
            print(f"file not found: {dest_file}")
        return dest_file

    def get_oss_cmd(self, subcommand="ls", options="", verbose=False):
        if self.is_valid_oss(verbose):
            # check ossutil command
            if not Path(self.ossutil).exists():
                msg = "\n".join(
                    [
                        f"ossutil command not found: {self.ossutil}",
                        f'find tosutil at: "https://www.volcengine.com/docs/6349/148777"',
                        f'find ossutil at: "https://help.aliyun.com/zh/oss/developer-reference/oss-tools"',
                    ]
                )
                print(msg)
                return None
            # tos require 'region'
            if self.oss_path.startswith("tos://"):
                opt_region = f"-re {self.region}"
            else:
                opt_region = ""
            # dest
            if subcommand in ["cp", "sync"]:
                opt_dest = self.out_dir
            else:
                opt_dest = ""
            # command
            cmd = " ".join(
                [
                    self.ossutil,
                    subcommand,
                    options,
                    opt_region,
                    f"-e {self.endpoint}",
                    f"-i {self.accesskeyid}",
                    f"-k {self.accesskeysecret}",
                    self.oss_path,
                    opt_dest,
                ]
            )
            return cmd

    def list_oss(self, verbose=False):
        cmd = self.get_oss_cmd("ls", verbose=verbose)
        if not isinstance(cmd, str):
            print(f"[failed] invalid oss_info")
            return 1  # exit code
        # 3. run
        result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
        # 4. check output
        if result.returncode > 0:
            print(result.stdout.decode())
            f_list = []
        else:
            s = result.stdout.decode("utf-8")
            # full size
            ## eg:
            ## 2023-10-18 14:01:50 +0800 CST 478748813 IA  D6BD...   oss://...
            ss = [
                re.search("\\s([0-9]+)\\s", i).group(1)
                for i in re.split("\\n", s)
                if bool(re.search("\\s([0-9]+)\\s", i))
            ]
            # file size
            f_size = sum(list(map(int, ss)))
            f_size = f"{f_size/1024**3:.1f}"  # GB
            # full_path of files, .gz
            f_list = [
                i.split()[-1] for i in re.split("\\n", s) if i.endswith("gz")
            ]
            msg_list = []
            if len(f_list) > 1:
                msg_files = [f_list[0], "...", f_list[-1]]
                msg_list = [
                    f"{i+1:>2}. {Path(k).name}"
                    for i, k in enumerate(msg_files)
                ]
            msg_list.append(
                f"[{len(f_list)}] fastq files, [{f_size} GB] in total"
            )
            msg = "\n".join(msg_list)
            print(msg)
        return [result.returncode, f_list]  # 0=yes, 1=no

    def file_exists(self):
        if self.is_valid_oss(verbose=False):
            _, f_list = self.list_oss()
            f_list = [Path(i).relative_to(self.oss_path) for i in f_list]
            ss = [i for i in f_list if Path(i).exists()]
            if len(ss) > 0:
                msg = "\\n".join(ss)
                msg += f"\\n[{len(ss)}] files found in out_dir: {self.out_dir}"
                print(msg)
            return len(ss) > 0

    def download_oss(self, subcommand="cp", verbose=False):
        """
        Download OSS files using ossutil/tosutil
        subcommand: cp, sync
        """
        # build command
        cmd = self.get_oss_cmd("cp", options="-r", verbose=verbose)
        cmd_txt = str(Path(self.out_dir) / "run.sh")
        try:
            with open(cmd_txt, "wt") as w:
                w.write(cmd + "\n")
        except Exception as e:
            print(f"failed, could not write to file: {e} {cmd_txt}")
            return None
        # show message
        msg = "\n".join(
            [
                "=" * 80,
                f'{"program":>14} : {"download_oss"}',
                f'{"date":>14} : {get_date()}',
                f'{"oss_path":>14}: {self.oss_path}',
                f'{"out_dir":>14}: {self.out_dir}',
                "=" * 80,
            ]
        )
        print(msg)
        # list files
        if self.file_exists() and not self.overwrite:
            print(f"file exists: {self.out_dir}")
        elif not self.dry_run:
            result = subprocess.run(cmd.split(), stderr=subprocess.PIPE)
            if result.returncode > 0:
                print(f"failed, ossutil {subcommand}, {result.stderr}")

    def run(self):
        if isinstance(self.http, str):
            self.download_http()
        else:
            self.download_oss()


class Download2(object):
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
            "out_dir": None,
            "http": None,
            "accesskeyid": None,
            "accesskeysecret": None,
            "oss_path": None,
            "endpoint": None,
            "region": None,
            "ossutil": None,
            "overwrite": False,
            "dry_run": False,
            "maxspeed": 10,  # KB/s
        }
        self = update_obj(self, args, force=False)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd() / "from_illumina")
        self.out_dir = str(Path(self.out_dir).expanduser().absolute())
        # self.is_http = isinstance(self.http, str)
        if isinstance(self.http, str) or self.is_valid_oss():
            if not Path(self.out_dir).exists():
                try:
                    Path(self.out_dir).mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(f"Could not create dir, {e}")

    def is_valid_oss(
        self,
    ):
        opts = {
            "accesskeyid": self.accesskeyid,
            "accesskeysecret": self.accesskeysecret,
            "oss_path": self.oss_path,
            "endpoint": self.endpoint,
            "region": self.region,
        }
        return all([isinstance(i, str) for i in list(opts.values())])

    def get_url_name(self, x):
        """
        get the filename from url
        http://xxx.com/kefu%2Fxxx2.tar?Expires=1xx&OSSAccessKeyId=xxx&Signature=xxx
        re.sub('^\w+\%|\?.*$', '', b)
        """
        if isinstance(x, str):
            xname = re.sub("?Expires=.*", "", Path(x).name)  # fix
            return re.sub("[^A-Za-z0-9\\_\\-\\.]", "", xname)
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
            print(f"skipped, http is not str, got {type(self.http)}")
            return None
        filename = self.get_url_name(self.http)
        dest_file = str(Path(self.out_dir) / filename)
        # show msg
        msg = "\n".join(
            [
                "=" * 80,
                f'{"Program":>14} : {"download_http"}',
                f'{"Date":>14} : {get_date()}',
                f'{"url":>14}: {self.download_http}',
                f'{"out_dir":>14}: {self.out_dir}',
                "=" * 80,
            ]
        )
        print(msg)
        # check_dir(self.out_dir)
        if Path(dest_file).exists() and not self.overwrite:
            print(f"file exists: {dest_file}")
        else:
            try:
                request.urlretrieve(self.http, dest_file)
            except Exception as e:
                print(f"failed downloading url, {e}, {self.http}")
        if not Path(dest_file).exists():
            print(f"file not found: {dest_file}")
        return dest_file

    def run(self):
        if self.is_valid_oss():
            args = self.__dict__.copy()  #
            if self.oss_path.startswith("tos"):
                TOS(**args).download()
            elif self.oss_path.startswith("oss"):
                OSS(**args).download()
            else:
                print(f"unknown oss_path: {self.oss_path}")
                pass
        elif isinstance(self.http, str):
            print("!A-3")
            self.download_http()
        else:
            pass


def download(**kwargs):
    """
    Download files using oss or http

    Args:
        oss_txt (str): A file contains oss info or http url
        out_dir (str): Ouput directory, default [pwd]

    Returns:
        None
    """
    for oss_info in OSSinfo(**kwargs).info:
        oss_info.update(
            {
                "out_dir": kwargs.get("out_dir", None),
                "overwrite": kwargs.get("overwrite", False),
                "dry_run": kwargs.get("dry_run", False),
                "maxspeed": kwargs.get("maxspeed", 10),  # 10MB/s
            }
        )
        Dnd = Download2
        ossutil = oss_info.get("ossutil", None)
        if isinstance(ossutil, str):
            if Path(ossutil).exists():
                Dnd = Download
        Dnd(**oss_info).run()


def get_args():
    example = "\n".join(
        [
            "Download fq data from OSS or http ...",
            "1. download oss data",
            "$ python download.py -s oss.txt -o out_dir",
        ]
    )
    parser = argparse.ArgumentParser(
        prog="download",
        description="download oss data",
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--oss-txt",
        dest="oss_txt",
        required=True,
        help="oss.txt file, could be a URL, or OSS config file",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        required=True,
        help="directory to save the results",
    )
    parser.add_argument(
        "-r",
        "--region",
        dest="region",
        default=None,
        help="choose region, eg: beijing, shanghai, ..., default [None]",
    )
    parser.add_argument(
        "-e",
        "--endpoint",
        dest="endpoint",
        default=None,
        help="specify the endpoint, eg: http://oss-cn-shanghai.aliyuncs.com, default [None]",
    )
    parser.add_argument(
        "-t",
        "--ossutil",
        dest="ossutil",
        default="/data/biosoft/ossutil/ossutil",
        help="ossutil command-line tool, default: [/data/biosoft/ossutil/ossutil]",
    )
    parser.add_argument(
        "-i",
        "--AccessKeyId",
        dest="accesskeyid",
        default=None,
        help="The AccessKeyId, default: [None]",
    )
    parser.add_argument(
        "-k",
        "--AccessKeySecret",
        dest="accesskeysecret",
        default=None,
        help="The AccessKeySecret, default: [None]",
    )
    parser.add_argument(
        "-O", "--overwrite", action="store_true", help="Overwrite exists files"
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Do not download files, just check the oss_info by ls command",
    )
    parser.add_argument(
        "--maxspeed",
        type=int,
        default=1000,
        help="set speedlimit KB/s, default [1000]",
    )
    return parser


def main():
    args = vars(get_args().parse_args())
    download(**args)


if __name__ == "__main__":
    main()
