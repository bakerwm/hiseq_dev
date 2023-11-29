#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Help functions for general purpose
# function
update_obj
run_shell_cmd
hash_string
check_hash_string
find_longest_common_str
gen_random_string
print_dict
init_cpu
get_date
is_cmd
unique_list
combination
convert_image
download_file
url_to_link

# class
Config
"""

__author__ = 'Ming Wang'
__email__ = 'wangm08 at hotmail.com'
__date__ = '2022-05-20'
# __version__ = '1.0.1'


import os
import sys
import json
import yaml
import toml
import pickle
import logging
import subprocess
import random
import string
import tempfile
import shutil
import hashlib # hash functions
# import uuid # generate a random number
from PIL import Image # pillow
from datetime import datetime
from dateutil import tz
from itertools import combinations
import Levenshtein as lev # distance
# import hiseq_dev
# from hiseq.utils.file import file_exists, list_dir
# from .helper import log
# from .file import file_exists, file_abspath
from difflib import SequenceMatcher
from urllib import request


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def update_obj(obj, d, force=True, remove=False):
    """
    Update the object, by dict
    Parameters:
    -----------
    obj : object
        an object
    d : dict
        a dict, save key:value pairs for object attribute:value
    force : bool
        update exists attributes to object, default: True
    remove : bool
        remove exists attributes from object, default: False
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


class Config(object):
    """
    Working with config, support formats: dict,yaml,yml,pickle
    Config().load(), and Config.dump()
    
    Example:
    1. write to file
    >>> Config(d).dump('out.json')
    >>> Config(d).dump('out.toml')
    >>> Config(d).dump('out.pickle')
    
    2. load from file
    >>> d = Config().load('in.yaml')
    """
    def __init__(self, x=None, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.x = x


    def load(self, x=None):
        """
        Read data from x, auto-recognize the file-type
        yaml
        toml
        json
        pickle
        txt
        ...
        """
        if x == None:
            x = self.x # dict or str
        if x is None:
            x_dict = None # {} ?
        elif isinstance(x, dict):
            x_dict = dict(sorted(x.items(), key=lambda i:i[0]))
        elif isinstance(x, str):
            reader = self.pick_reader(x)
            if reader is None:
                x_dict = None
                log.error('unknown x, {}'.format(x))
            else:
                x_dict = reader(x)
        else:
            x_dict = None
            log.warning('dump(x=) dict,str expect, got {}'.format(
                type(x).__name__))
        return x_dict


    def dump(self, d=None, x=None):
        if d is None:
            d = self.load(self.x)
        # make sure: dict
        if isinstance(x, str):
            writer = self.pick_writer(x)
            if writer is None:
                log.error('unknown x, {}'.format(x))
            else:
                writer(d, x)
        else:
            log.warning('dump(x=) expect str, got {}'.format(
                type(x).__name__))


    def guess_format(self, x):
        """
        Guess the file format by file extension    
        str: yaml, yml, toml, json, pickle
        dict: dict
        """
        fmts = {
            'json': 'json',
            'yaml': 'yaml',
            'yml': "yaml",
            'toml': 'toml',
            'pickle': 'pickle',
        }
        if isinstance(x, dict):
            fmt = 'dict'
        elif isinstance(x, str):
            ext = os.path.splitext(x)[1]
            ext = ext.lstrip('.').lower()
            fmt = fmts.get(ext, None)
        else:
            fmt
        return fmt


    def pick_reader(self, x):
        fmt = self.guess_format(x)
        readers = {
            'json': self.from_json,
            'yaml': self.from_yaml,
            'toml': self.from_toml,
            'pickle': self.from_pickle
        }
        return readers.get(fmt, None)


    def pick_writer(self, x):
        fmt = self.guess_format(x)
        writers = {
            'json': self.to_json,
            'yaml': self.to_yaml,
            'toml': self.to_toml,
            'pickle': self.to_pickle
        }
        return writers.get(fmt, None)

    
    # json 
    def from_json(self, x):
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = json.load(r) # sorted by key
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_json() failed, {}'.format(exc))
        else:
            log.error('from_json() failed, file not exists: {}'.format(x))
        return d


    def to_json(self, d, x):
        """
        Save dict to file as json format
        """
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_json(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_json(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_json(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    json.dump(d, w, indent=4, sort_keys=True)
                # return x
            except Exception as exc:
                log.error('to_json() failed, {}'.format(exc))


    # YAML
    def from_yaml(self, x):
        d = None
        if self.guess_format(x) == 'yaml':
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = yaml.load(r, Loader=yaml.FullLoader)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_yaml() failed, {}'.format(exc))
        else:
            log.error('from_yaml() failed, not a yaml file: {}'.format(x))
        return d


    def to_yaml(self, d, x):
        """
        Save dict to file as YAML format
        yaml.dump(), does not support OrderedDict
        Solution: OrderedDict -> json -> dict
        """
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_yaml(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_yaml(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_yaml(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    yaml.dump(d, w)
            except:
                log.warning('saving as YAML failed, use TOML instead')
                x_toml = os.path.splitext(x)[0] + '.toml'
                with open(x_toml, 'wt') as w:
                    toml.dump(d, w)                
#             except Exception as exc:
#                 log.error('to_yaml() failed, {}'.format(exc))

    # TOML
    def from_toml(self, x):
        d = None
        if self.guess_format(x) == 'yaml':
            try:
                with open(x, 'r') as r:
                    if os.path.getsize(x) > 0:
                        d = toml.load(x)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_toml() failed, {}'.format(exc))
        else:
            log.error('from_toml() failed, file not exists: {}'.format(x))
        return d


    def to_toml(self, d, x):
        """
        Save dict to file as TOML format
        """        
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_toml(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_toml(d=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_toml(d=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wt') as w:
                    toml.dump(d, w)
                # return x
            except Exception as exc:
                log.error('to_toml() failed, {}'.format(exc))

    # pickle
    def from_pickle(self, x):
        d = None
        if os.path.exists(x):
            try:
                with open(x, 'rb') as r:
                    if os.path.getsize(x) > 0:
                        d = pickle.load(r)
                        d = dict(sorted(d.items(), key=lambda x:x[0]))
                        # d = collections.OrderedDict(sorted(d.items()))
            except Exception as exc:
                log.error('from_pickle() failed, {}'.format(exc))
            finally:
                return d
        else:
            log.error('from_pickle() failed, file not exists: {}'.format(x))

    def to_pickle(self, d, x):
        """
        Writing data to pickle file
        d dict, data to file
        x str, path to pickle file
        """        
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_pickle(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_pickle(x=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_pickle(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                with open(x, 'wb') as w:
                    pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)
                # return x
            except Exception as exc:
                log.error('to_pickle() failed, {}'.format(exc))

    
    # log, plain text file
    def to_log(self, d, x, stdout=False):
        """
        Writing data to log file: key: value format
        d dict, data to file
        x str, path to pickle file
        """
        x = os.path.abspath(x)
        if not isinstance(d, dict):
            log.error('to_log(d=) failed, dict expect, got {}'.format(
                type(d).__name__))
        elif not isinstance(x, str):
            log.error('to_log(x=) failed, str expect, got {}'.format(
                type(x).__name__))
        elif not os.path.exists(os.path.dirname(x)):
            log.error('to_log(x=) failed, file not exists: {}'.format(x))
        else:
            try:
                # organize msg
                msg = []
                for k, v in d.items():
                    if isinstance(v, str) or isinstance(v, numbers.Number) or isinstance(v, bool):
                        v = str(v)
                    elif isinstance(v, list):
                        v = ', '.join(map(str, v))
                    else:
                        v = '...' # skip
                    msg.append('{:30s} | {:<40s}'.format(k, v))
                # save
                with open(x, 'wt') as w:
                    w.write('\n'.join(msg) + '\n')
                if stdout:
                    print('\n'.join(msg))
                # return x
            except Exception as exc:
                log.error('to_log() failed, {}'.format(exc))


    def _tmp(self, suffix='.txt'):
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
            delete=False)
        return tmp.name


def run_shell_cmd(cmd):
    """
    This command is from 'ENCODE-DCC/atac-seq-pipeline'
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_common.py
    save log to file
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    cmd_name = '{} ...'.format(os.path.basename(cmd.split()[0]))
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd_name))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid,
        pgid,
        rc,
        stderr.strip(),
        stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            log.error(err_str)
    return (rc, stdout.strip('\n'), stderr.strip('\n'))


################################################################################
# ## for hash string; version with random number
# see: https://www.pythoncentral.io/hashing-strings-with-python/
#
# equivalent in R
# > digest::digest(s, algo = "sha256", serialize=FALSE, raw=FALSE)

# def hash_string(s):
#     salt = uuid.uuid4().hex # random number
#     return hashlib.sha256(salt.encode() + s.encode()).hexdigest() + ':' + salt

# def check_hashed_string(hash_s, s):
#     hs, salt = hash_s.split(':')
#     return hs == hashlib.sha256(salt.encode() + s.encode()).hexdigest()
 

## for hash string; version with random number
def hash_string(s):
    return hashlib.sha256(s.encode()).hexdigest()


def check_hash_string(hash_s, s):
    """
    Parameters:
    hash_s  : str
        The SHA-256 value for the input string 
        Also, could be the first few characters of the SHA-256 value
        Suggest useing the full version (64 characters) for checking
        
    s  : str
        The string for checking
    """
    # hash_s, could be the first few characters, (8?)
    hs = hash_string(s)
    return hs.startswith(hash_s) or hs == s


################################################################################
## tmp functions
def str_distance(x, y, partial=True):
    if isinstance(x, str) and isinstance(y, str):
        try:
            if partial:
                x = x[:len(y)]
                y = y[:len(x)]
            out = lev.distance(x, y)
        except:
            out = -1 # huge number
    else:
        out = -1
    return out


def find_longest_common_str(s1, s2):
    if isinstance(s1, str) and isinstance(s2, str):
        m = SequenceMatcher(None, s1, s2) # match
        l = m.find_longest_match(0, len(s1), 0, len(s2))
        out = s1[l.a:(l.a+l.size)]
    else:
        log.error('only support str, got s1={} s2={}'.type(
            type(s1).__name__, type(s2).__name__))
        out = None
    return out


def gen_random_string(slen=10):
    return ''.join(random.sample(string.ascii_letters + string.digits, slen))


def print_dict(d):
    d = dict(sorted(d.items(), key=lambda x:x[0]))
    # d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def init_cpu(threads=1, parallel_jobs=1):
    """
    The number of threads, parallel_jobs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()
    max_jobs = int(n_cpu / 4.0)
    ## check parallel_jobs (max: 1/4 of n_cpus)
    if parallel_jobs > max_jobs: 
        log.warning('Too large, change parallel_jobs from {} to {}'.format(
            parallel_jobs, max_jobs))
        parallel_jobs = max_jobs
    ## check threads
    max_threads = int(0.8 * n_cpu / parallel_jobs)
    if threads * parallel_jobs > 0.8 * n_cpu:
        log.warning('Too large, change threads from {} to {}'.format(
            threads, max_threads))
        threads = max_threads
    return (threads, parallel_jobs)


def get_date(timestamp=False):
    """
    Return the current date in UTC.timestamp or local-formated-string
    calculation in UTC
    return in local
    
    switch between timezone by datautil.tz
    
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
    
    Arguments
    ---------
    ts:  float
        The timestamp (microseconds), if not, get the current time
    
    in_seconds:  bool
        Return the date in seconds (microseconds?!)
    """
    now = datetime.now(tz.tzlocal())
    if isinstance(timestamp, bool) and timestamp:
        out = now.timestamp()
    else:
        now = now.astimezone(tz.tzlocal()) # to local
        out = now.strftime('%Y-%m-%d %H:%M:%S') # YY-mm-dd H:M:S
    return out


def is_cmd(x):
    """Check if the executable command"""
    return lambda i: shutil.which(i) is not None


def which(command):
    """
    shutil.which()
    subprocess.call()
    see: https://stackoverflow.com/a/377028/2530783
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(command)
    if fpath:
        if is_exe(command):
            return command
    else:
        for path in os.environ.get("PATH", "").split(os.pathsep):
            exe_file = os.path.join(path, command)
            if is_exe(exe_file):
                return exe_file

    return None


def unique_list(seq, sorted=True, idfun=None):
    """
    Get the unique items of list
    seq: a list with items
    sorted: whether sort the output(unique)
    get the unique of inlist
    see1: Markus
    remove duplicates from a list while perserving order
    https://stackoverflow.com/a/480227/2530783
    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    see2: ctcherry
    https://stackoverflow.com/a/89202/2530783
    def f5(seq, idfun=None):
        # order preserving
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            # in old Python versions:
            # if seen.has_key(marker)
            # but in new ones:
            if marker in seen: continue
            seen[marker] = 1
            result.append(item)
        return result
    """
    if idfun is None:
        def idfun(x): return x # for None
    if not isinstance(seq, list):
        log.error('list required, but get {}'.format(type(seq)))
        out = [] # blank
    elif sorted is True:
        out = list(set(seq))
    else:
        seen = set()
        out = [x for x in seq if x not in seen and not seen.add(x)]
    return out


def combination(x, n=2, return_index=True):
    """
    Generate the combination for a list of items    
    Parameters
    ----------
    x : list 
        A list of items
    example
    >>> x = ['ctl', 'ctl', 'exp', 'exp'] # keep order
    >>> combination(x, return_index=True)
    [[0, 1], [2, 3]]
    >>> combination(x, return_index=False)
    ['ctl', 'exp']
    """
    xu = unique_list(x, sorted=False) # keep order
    out = []
    if len(xu) >= n:
        item_pairs = list(combinations(xu, n))
        # for index
        index_pairs = []
        for (a, b) in item_pairs:
            index_a = [i for i, j in enumerate(x) if j == a]
            index_b = [i for i, j in enumerate(x) if j == b]
            index_pairs.append([index_a, index_b])
        out = index_pairs if return_index else item_pairs
    return out


def convert_image(x, out_fmt='PNG'):
    if not out_fmt in ['PNG', 'JPEG', "TIFF"]:
        log.error('out_fmt: [PNG|JPEG|TIFF], {} got'.format(out_fmt))
    out_ext = out_fmt.lower()
    out_img = os.path.splitext(x)[0] + '.' + out_ext
    # read/write
    if os.path.exists(out_img):
        log.warning('file exists, skipping ...: {}'.format(out_img))
    else:
        img = Image.open(x)
        img.save(out_img, out_fmt)

        
def download_file(url, file):
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
    file_dir = os.path.dirname(file)
    if not os.path.exists(file_dir):
        log.error('dir not exists: {}'.format(file_dir))
    elif os.path.exists(file):
        log.error('file exists: {}'.format(file))
    else:
        try:
            request.urlretrieve(url, file)
        except:
            log.error('failed downloading file: {}'.format(url))
        if os.path.exists(file):
            log.info('saving file: {}'.format(file))

        
def url_to_link(url, name=None, format='markdown'):
    """
    Convert url to link, in different format
    markdown, html
    
    markdown: [name](url) 
    html: <a href=url target='_blank'>name</a>
    """
    # support str, list
    if isinstance(url, str):
        if not isinstance(name, str):
            name = url
        if format == 'markdown':
            out = '[{name}]({url})'.format(name=name, url=url)
        elif format == 'html':
            out = "<a href={url} target='_blank'>{name}</a>".format(name=name, url=url)
        else:
            out = url
    elif isinstance(url, list):
        if not (isinstance(name, list) and len(url) == len(name)):
            name = url
        out = [url_to_link(i, n, format) for i,n in zip(url, name)]
    else:
        log.error('url illegal, str, list expect, got {}'.format(type(url).__name__))
        out = None
    return out


## load package files
## deprecated: failed, could not import 
# def list_pkg_file(*args):
#     """
#     list the file in package
#     >>> list_pkg_file('hiseq.data', 'illumina_index.csv')
#     see importlib.resources.path()
#     """
#     try:
#         with res.path(*args) as f:
#             out = str(f)
#     except:
#         print('could not find file: {}'.format(args[-1]))
#         out = None
#     return out
