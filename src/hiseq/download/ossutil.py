#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download files from different sources:
1. OSS (aliyun, ossutil)
2. TOS (volcengine, tosutil)
2. url
"""

# import os
# import re
# import argparse
# import subprocess
import sys
from pathlib import Path
from urllib import request
from dateutil import tz
from datetime import datetime

# from hiseq.utils.utils import update_obj, get_date

import time
import tos
from tos import DataTransferType, DownloadEventType
from tos.checkpoint import CancelHook

import oss2
from oss2.models import OSS_TRAFFIC_LIMIT


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


class TOS(object):
    """
    Access tos objects using python SDK `tos`
    see documentation at: https://www.volcengine.com/docs/6349/92785

    Examples:
    >>> ss = TOS(ak, sk, oss_path, endpoint, region)
    >>> ss.ls()
    >>> ss.download()
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args = {
            "accesskeyid": None,
            "accesskeysecret": None,
            "oss_path": None,
            "endpoint": None,
            "region": None,
            "out_dir": None,
            "bucket_name": None,  # parse from oss_path
            "prefix": None,  # parse from oss_path
            "overwrite": False,
            "dry_run": False,
            "maxspeed": 1000,  # 1000 KB/s, ~1 MB/s
        }
        self = update_obj(self, args, force=False)
        parts = self.parse_oss_path()
        if isinstance(parts, list):
            self.oss_path, self.bucket_name, self.prefix = parts[:3]
        self.client = self.open_client()
        self.is_valid_oss = isinstance(self.client, tos.TosClientV2)
        # output dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd() / "from_illumina")
        self.out_dir = str(Path(self.out_dir).expanduser().absolute())
        if self.is_valid_oss:
            if not Path(self.out_dir).exists():
                try:
                    Path(self.out_dir).mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(f"Could not create dir, {e}")

    def parse_oss_path(self, x=None):
        """
        Parse bucket_name, prefix from oss_path
        oss_path: 'oss://seekgene-release/2307015_1010/release/'
        """
        if x is None:
            x = self.oss_path
        if isinstance(x, str):
            if x[2:].startswith("s://"):
                # oss_path =  x.strip('/') + '/' # update x, add suffix
                oss_path = x
                parts = Path(oss_path).parts
                bucket_name = parts[1]  # seekgene-release
                url_base = f"{parts[0]}//{parts[1]}"
                prefix = str(Path(oss_path).relative_to(url_base))
                return [oss_path, bucket_name, prefix]

    def open_client(self):
        try:
            client = tos.TosClientV2(
                self.accesskeyid,
                self.accesskeysecret,
                self.endpoint,
                self.region,
            )
            return client
        except tos.exceptions.TosClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e.message}, cause: {e.cause}")
        except tos.exceptions.TosServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with message: {e.message}")
            print(f"error with http code: {e.status_code}")
            print(f"error with ec: {e.ec}")
            print(f"error with request url: {e.request_url}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def close_client(self):
        if self.is_valid_oss:
            self.client.close()

    def ls_content(self, prefix=None):
        if prefix is None:
            prefix = self.prefix  #
        if not self.is_valid_oss:
            return None
        try:
            content_list = []
            is_truncated = True
            next_continuation_token = ""
            while is_truncated:
                out = self.client.list_objects_type2(
                    self.bucket_name,
                    prefix=prefix,
                    delimiter="/",
                    continuation_token=next_continuation_token,
                )
                is_truncated = out.is_truncated
                next_continuation_token = out.next_continuation_token

                # recursive list sub-folders
                for px in out.common_prefixes:
                    # print(f'subdir: {px.prefix}')
                    content_list.extend(self.ls_content(px.prefix))
                    # break #!!! debug
                # list files under prefix
                for content in out.contents:
                    # print(f'file: {content.key}')
                    content_list.append(content)
                    # break #!!! debug
            return content_list
        except tos.exceptions.TosClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e.message}, cause: {e.cause}")
        except tos.exceptions.TosServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with message: {e.message}")
            print(f"error with http code: {e.status_code}")
            print(f"error with ec: {e.ec}")
            print(f"error with request url: {e.request_url}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def ls(self):
        content_list = self.ls_content()
        if isinstance(content_list, list):
            key_list = [content.key for content in content_list]
            size_list = [content.size for content in content_list]
            # message
            total_size = sum(size_list) / 1024**3  # GB
            if len(key_list) > 0:
                print(f"{key_list[0]}\n...")
                print(
                    f"[{len(key_list)}] files [{total_size:.1f}] GB in total"
                )
            return key_list

    def download_key(self, key):
        try:
            if isinstance(key, str):
                #!!! caution, key contains prefix only
                key_path = f"tos://{self.bucket_name}/{key}"
                local_file = str(
                    Path(self.out_dir)
                    / Path(key_path).relative_to(self.oss_path)
                )
                if Path(local_file).exists() and not self.overwrite:
                    print(f"file exists: {local_file}")
                else:
                    self.client.download_file(
                        self.bucket_name,
                        key,
                        local_file,
                        task_num=3,
                        part_size=1 * 1024 * 1024,
                        # download_event_listener=self.utils('download_event'),
                        data_transfer_listener=self.utils("percentage"),
                        rate_limiter=self.utils("rate_limiter"),
                    )
                return local_file
            else:
                print(f"download(key) expect str, got {type(key)}")
                return None
        except tos.exceptions.TosClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e.message}, cause: {e.cause}")
        except tos.exceptions.TosServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with message: {e.message}")
            print(f"error with http code: {e.status_code}")
            print(f"error with ec: {e.ec}")
            print(f"error with request url: {e.request_url}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def download(self):
        key_list = self.ls()
        if isinstance(key_list, list):
            if not bool(self.dry_run):
                [self.download_key(key) for key in key_list]

    def upload(self):
        pass

    def utils(self, name="percentage"):
        # progress bar
        def percentage(
            consumed_bytes, total_bytes, rw_once_bytes, type: DataTransferType
        ):
            if total_bytes:
                rate = float(consumed_bytes) / float(total_bytes) * 100
                msg = ", ".join(
                    [
                        f"rate:{rate:.0f}%",
                        f"consumed_bytes:{consumed_bytes}",
                        f"total_bytes:{total_bytes}",
                        f"rw_once_bytes:{rw_once_bytes}",
                        f"type:{type}",
                    ]
                )
                print(msg, end="\r")

        # event
        def download_event(
            type: DownloadEventType,
            err,
            bucket,
            key,
            version_id,
            file_path,
            checkpint_file,
            tmp_file,
            download_part,
        ):
            print(
                type,
                err,
                bucket,
                key,
                version_id,
                file_path,
                checkpint_file,
                tmp_file,
                download_part,
                end="\r",
            )

        # set speed limitation, eg: 10 KB/s, max 10 + 5 KB/s
        rate_limiter = tos.RateLimiter(
            rate=self.maxspeed * 1024 * 8, capacity=5 * 1024 * 8
        )

        # # 继承 CancelHook 类实现断点续传下载任务取消功能
        # class MyCancel(CancelHook):
        #     def cancel(self, is_abort: bool):
        #         # is_abort 为 true 时删除上下文信息并 abort 分段上传任务，
        #         # 为 false 时只是中断当前执行
        #         # 重写 cancel 方法时必须调用 父类的 cancel 方法
        #         # 模拟 10 秒后取消任务
        #         time.sleep(10)
        #         super(MyCancel, self).cancel(is_abort=is_abort)
        #         print('some user define')
        # cancel = MyCancel()

        # output
        units = {
            "percentage": percentage,
            "download_event": download_event,
            "rate_limiter": rate_limiter,
            # 'cancel': MyCancel,
        }

        return units.get(name, None)


class OSS(object):
    """
    Access tos objects using python SDK `oss2`
    see documentation at: https://help.aliyun.com/zh/oss/developer-reference/python

    Examples:
    >>> ss = OSS(ak, sk, oss_path, endpoint, region)
    >>> ss.ls()
    >>> ss.download()
    """

    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args = {
            "accesskeyid": None,
            "accesskeysecret": None,
            "oss_path": None,
            "endpoint": None,
            "region": None,
            "out_dir": None,
            "bucket_name": None,  # parse from oss_path
            "prefix": None,  # parse from oss_path
            "overwrite": False,
            "dry_run": False,
            "maxspeed": 1000,  # 1000 KB/s, ~1 MB/s
        }
        self = update_obj(self, args, force=False)
        parts = self.parse_oss_path()
        if isinstance(parts, list):
            self.oss_path, self.bucket_name, self.prefix = parts[:3]
        self.bucket = self.open_bucket()
        self.is_valid_oss = isinstance(self.bucket, oss2.Bucket)
        # output dir
        if not isinstance(self.out_dir, str):
            self.out_dir = str(Path.cwd() / "from_illumina")
        self.out_dir = str(Path(self.out_dir).expanduser().absolute())
        if self.is_valid_oss:
            if not Path(self.out_dir).exists():
                try:
                    Path(self.out_dir).mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    print(f"Could not create dir, {e}")

    def parse_oss_path(self, x=None):
        """
        Parse bucket_name, prefix from oss_path
        oss_path: 'oss://seekgene-release/2307015_1010/release/'
        """
        if x is None:
            x = self.oss_path
        if isinstance(x, str):
            if x[2:].startswith("s://"):
                # oss_path =  x.strip('/') + '/' # update x, add suffix
                oss_path = x
                parts = Path(oss_path).parts
                bucket_name = parts[1]  # seekgene-release
                url_base = f"{parts[0]}//{parts[1]}"
                prefix = str(Path(oss_path).relative_to(url_base))
                prefix = prefix.rstrip("/") + "/"
                return [oss_path, bucket_name, prefix]

    def open_bucket(self):
        try:
            bucket = oss2.Bucket(
                oss2.Auth(self.accesskeyid, self.accesskeysecret),
                self.endpoint,
                self.bucket_name,
            )
            return bucket
        except oss2.exceptions.ClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e}")
        except oss2.exceptions.ServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with ec: {e.ec}")
            print(f"error with http code: {e.status}")
            print(f"error with message: {e.message}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def ls_content(self, prefix=None):
        if prefix is None:
            prefix = self.prefix  #
        if not self.is_valid_oss:
            return None
        try:
            content_list = []
            for obj in oss2.ObjectIterator(
                self.bucket, prefix=prefix, delimiter="/"
            ):
                if obj.is_prefix():
                    # sub-folders
                    # print(f'sub_dir: {obj.key}')
                    content_list.extend(self.ls_content(obj.key))
                else:
                    # files
                    ## !!! in-case 'path/'
                    if obj.key.endswith("/"):
                        continue
                    # print(f'file: {obj.key}')
                    content_list.append(obj)
            return content_list
        except oss2.exceptions.ClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e}")
        except oss2.exceptions.ServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with ec: {e.ec}")
            print(f"error with http code: {e.status}")
            print(f"error with message: {e.message}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def ls(self):
        content_list = self.ls_content()
        if isinstance(content_list, list):
            key_list = [content.key for content in content_list]
            size_list = [content.size for content in content_list]
            # message
            total_size = sum(size_list) / 1024**3  # GB
            if len(key_list) > 0:
                print(f"{key_list[0]}\n...")
                print(
                    f"[{len(key_list)}] files [{total_size:.1f}] GB in total"
                )
            return key_list

    def download_key(self, key):
        try:
            if isinstance(key, str):
                #!!! caution, key contains prefix only
                key_path = f"oss://{self.bucket_name}/{key}"
                local_file = str(
                    Path(self.out_dir)
                    / Path(key_path).relative_to(self.oss_path)
                )
                if Path(local_file).exists() and not self.overwrite:
                    # print(f'file exists: {local_file}')
                    pass
                else:
                    # local directory
                    local_dir = Path(local_file).parent
                    if not Path(local_dir).exists():
                        Path(local_dir).mkdir(parents=True, exist_ok=True)
                    # speedlimit, maxspeed
                    maxspeed = self.maxspeed * 1024 * 8  # KB/s
                    headers = dict()
                    headers[OSS_TRAFFIC_LIMIT] = str(maxspeed)
                    # run
                    oss2.resumable_download(
                        self.bucket,
                        key,
                        local_file,
                        store=oss2.ResumableDownloadStore(root="/tmp"),
                        multiget_threshold=1
                        * 1024
                        * 1024,  # threshold for multiget, 1 MB
                        part_size=1 * 1024 * 1024,  # 1 MB [100 KB - 5 GB]
                        progress_callback=self.utils("percentage"),
                        num_threads=3,
                        headers=headers,
                    )
                return local_file
            else:
                print(f"download(key) expect str, got {type(key)}")
                return None
        except oss2.exceptions.ClientError as e:
            # catch error from client, illegal args
            print(f"fail with client, message:{e}")
        except oss2.exceptions.ServerError as e:
            print(f"fail with server error, code: {e.code}")
            print(f"error with request id: {e.request_id}")
            print(f"error with ec: {e.ec}")
            print(f"error with http code: {e.status}")
            print(f"error with message: {e.message}")
        except Exception as e:
            print(f"fail with unknown error: {e}")

    def download(self):
        key_list = self.ls()
        if isinstance(key_list, list):
            if not bool(self.dry_run):
                [self.download_key(key) for key in key_list]
        sys.stdout.flush()  # clear stdout, screen

    def upload(self):
        pass

    def utils(self, name="percentage"):
        # progress bar
        def percentage(consumed_bytes, total_bytes):
            if total_bytes:
                rate = float(consumed_bytes) / float(total_bytes) * 100
                msg = ", ".join(
                    [
                        f"rate:{rate:.0f}%",
                        f"consumed_bytes:{consumed_bytes}",
                        f"total_bytes:{total_bytes}",
                        f"type:{type}",
                    ]
                )
                print(msg, end="\r")
                # print('\r{0}% '.format(rate), end='')
                sys.stdout.flush()

        # output
        units = {"percentage": percentage}

        return units.get(name, None)
