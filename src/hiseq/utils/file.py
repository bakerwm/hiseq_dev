#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
functions for file manipulation
## 1. manipulate file
- check_dir
- copy_dir
- remove_dir
- check_file
- copy_file
- remove_file
- rename_file
- symlink_file
- gzip_file
- read_lines => read_file ?
## 2. file info
- file_is_gzipped
- file_prefix
- file_abspath
- file_exists
- file_nrows
## 3. search files
- list_dir
- list_file
"""

import os
#  import sys
# import re
import pathlib
import shutil
# import logging
import fnmatch
import binascii
# import subprocess
# import pyfastx
import pysam
import pyBigWig
from xopen import xopen
# from hiseq.utils.seq import Fastx
from hiseq.utils.utils import log


## 1. manipulate files
def check_dir(x, **kwargs):
    """Check if x is path
    Parameters
    ----------
    x : str
        Path to a file
    Keyword Parameters
    ------------------
    show_error : bool
        Show the error messages
    show_log : bool
        Show the log messages
    create_dirs : bool
        Create the dirs
    """
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    create_dirs = kwargs.get('create_dirs', True) # default: True
    if isinstance(x, str):
        out = False
        if os.path.isdir(x):
            out = True
        elif os.path.isfile(x):
            if show_error:
                log.error('file exists, not a directory: {}'.format(x))
        else:
            if create_dirs:
                try:
                    os.makedirs(x)
                    out = True
                except:
                    if show_error:
                        log.error('`os.makedirs` failed: {}'.format(x))
        # show log
        flag = 'ok' if out else 'failed'
        if show_log is True:
            log.info('{:<6s} : {}'.format(flag, x))
    elif isinstance(x, list):
        out = all([check_dir(i, **kwargs) for i in x])
    else:
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
    return out


def copy_dir(src, dest, force=False):
    """
    Copy the whole directory
    Parameters
    ----------
    x : str or list
        The file(s) to be removed
    force : bool
        Copy files, overwrite dest file
    """
    if isinstance(src, str) and isinstance(dest, str):
        if os.path.isdir(src):
            shutil.copytree(src, dest)
        else:
            log.error('src is not directory')
    else:
        log.error('both src and dest required str, got {}, {}'.format(
            type(src).__name__, type(dest).__name__))


def remove_dir(x, **kwargs):
    """
    Remove directory
    Parameters
    ----------
    x : str or list
        The file(s) to be removed
    ask : bool
        Ask the user, before the files removed.
    check_empty : True
        Do not proceed, if directory is not empty
    show_log : True
        Display the status of the files
    show_error : False
        Display the error messages
    """
    ask = kwargs.get('ask', True)
    check_empty = kwargs.get('check_empty', True)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    if isinstance(x, str):
        if os.path.isdir(x):
            x_files = list_dir(x, full_names=True, recursive=False, include_dirs=True)
            is_empty = len(x_files) == 0
            if check_empty and not is_empty:
                is_rm = False
            else:
                is_rm = True
            # rm, ask
            if is_rm:
                ask_msg = input('Remove: {}, [Y|n]: '.format(x)) if ask else 'Y'
            else:
                ask_msg = 'no'
            rm_tag = 'yes' if ask_msg.lower() in ['y', 'yes'] else 'no'
            empty_tag = 'yes' if is_empty else 'no'
            # do-the-thing, removing
            if rm_tag == 'yes':
                try:
                    shutil.rmtree(x)
                except:
                    rm_tag = 'no'
                    if show_error:
                        log.error('failed, remove path: {}'.format(x))
        else:
            rm_tag = 'no'
            empty_tag = 'NA'
            if show_error:
                log.error('x is not path, {}'.format(x))
        if show_log:
            log.info('rm:{:3s}\tis_empty:{}\t{:3s}'.format(rm_tag, empty_tag, x))
    elif isinstance(x, list):
        [remove_dir(i, **kwargs) for i in x]
    elif isinstance(x, dict):
        for k,v in x.items:
            if isinstance(v, str) or isinstance(v, list):
                remove_dir(v, **kwargs)
    else:
        log.error('x, str or list or dict expected, got {}'.format(
            type(x).__name__))


def check_file(x, **kwargs):
    """
    Check the x file
    1. file exists
    Parameters
    ----------
    x : str
        Path to a file
    Keyword Parameters
    ------------------
    show_error : bool
        Show the error messages
    show_log : bool
        Show the log messages
    check_empty : bool
        Check if the file is empty or not,  gzipped empty file, size=20
    emptycheck : bool
        see check_empty
    """
    show_error = kwargs.get('show_error', False)
    show_log = kwargs.get('show_log', False)
    check_empty = kwargs.get('check_empty', False)
    if isinstance(x, str):
        if file_exists(x):
            x_size = os.stat(x).st_size
            # empty gzipped file, size=20
            q_size = 20 if x.endswith('.gz') else 0
            out = x_size > q_size if check_empty else True
            if show_log:
                flag = 'ok' if out else 'failed'
                log.info('{:<6s} : {}'.format(flag, x))
        else:
            if show_error:
                log.error('file not exists: {}'.format(x))
            out = False # failed
    elif isinstance(x, list):
        out = all([check_file(i, **kwargs) for i in x])
    else:
        if show_error:
            log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out


def copy_file(src, dest, force=False):
    """Copy file

    Parameters
    ----------
    x : str or list
        The file(s) to be removed

    force : bool
        Copy files, overwrite dest file
    """
    if not isinstance(src, str):
        log.error('src, expect str, got {}'.format(type(src).__name__))
    elif not isinstance(dest, str):
        log.error('dest, expect str, got {}'.format(type(src).__name__))
    elif os.path.isfile(src):
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
        else:
            dest_file = dest
        # do-the-thing
        if file_exists(dest_file) and not force:
            log.error('copy_file() skipped, dest exists: {}'.format(dest_file))
        else:
            try:
                shutil.copy(src, dest_file)
            except:
                log.error('copy_file() failed, {}'.format(dest_file))
    else:
        log.warning('copy_file() failed, src not vaild: {}'.format(src))


def remove_file(x, **kwargs):
    """Remove files

    Parameters
    ----------
    x : str or list
        The file(s) to be removed

    ask : bool
        Ask the user, before the files removed.
    """
    ask = kwargs.get('ask', True)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    if isinstance(x, str):
        if os.path.isfile(x):
            file_tag = 'yes'
            ask_msg = input('Remove: {}, [Y|n]: '.format(x)) if ask else 'Y'
        else:
            file_tag = 'no'
            ask_msg = 'no'
        rm_tag = 'yes' if ask_msg.lower() in ['y', 'yes'] else 'no'
        # do-the-thing, removing
        if rm_tag == 'yes':
            try:
                os.remove(x)
            except:
                rm_tag = 'no'
                if show_error:
                    log.error('failed, remove file: {}'.format(x))
        if show_log:
            log.info('rm:{:3s}\tis_file:{}\t{:3s}'.format(rm_tag, file_tag, x))
    elif isinstance(x, list):
        [remove_file(i, **kwargs) for i in x]
    elif isinstance(x, dict):
        for k,v in x.items:
            if isinstance(v, str) or isinstance(v, list):
                remove_file(v, **kwargs)
    else:
        log.error('x, str or list or dict expected, got {}'.format(
            type(x).__name__))


def symlink_file(src, dest, absolute_path=False, force=False):
    """
    Create symlink
    Parameters
    ----------
    src : str
        The source file
    dest : str
        The target file, dir exists
    absolute_path : bool
        Use abs_path instead
    force : bool
        Copy files, overwrite dest file
    """
    if not isinstance(src, str):
        log.error('src, expect str, got {}'.format(type(src).__name__))
    elif not isinstance(dest, str):
        log.error('dest, expect str, got {}'.format(type(src).__name__))
    elif os.path.isfile(src):
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
        else:
            dest_file = dest
        dest_file = os.path.abspath(os.path.expanduser(os.path.expandvars(dest_file)))
        # the relative path of src
        src_dir = os.path.dirname(src)
        src_name = os.path.basename(src)
        dest_dir = os.path.dirname(dest_file)
        src_dir_rel = os.path.relpath(src_dir, dest_dir)
        src_rel = os.path.join(src_dir_rel, src_name)
        src_file = src if absolute_path else src_rel
        # do-the-thing
        if file_exists(dest_file) and not force:
            log.info('symlink_file() skipped, dest exists: {}'.format(dest_file))
        else:
            try:
                os.symlink(src_file, dest_file)
            except:
                log.error('symlink_file() failed, {}'.format(dest_file))
    elif os.path.islink(src):
        pass
    else:
        log.warning('symlink_file() failed, src not vaild: {}'.format(src))


def gzip_file(src, dest=None, decompress=True, **kwargs):
    """Gzip Compress or Decompress files using gzip module in python
    rm, True/False, whether remove old file

    # check the src file by extension: .gz
    """
    a = os.path.exists(src)
    b = file_exists(src)
    show_log = kwargs.get('show_log', True)
    show_error = kwargs.get('show_error', False)
    compresslevel = kwargs.get('compresslevel', 1)
    threads = kwargs.get('threads', 4)
    flag = False
    # input: src
    if not isinstance(src, str):
        if show_error:
            log.error('src expect str, got {}'.format(src))
    if not file_exists(src):
        if show_error:
            log.error('src not exists, {}'.format(src))
    # output: dest
    if dest is None:
        dest = os.path.splitext(src)[0] if decompress else src + '.gz'
    if isinstance(dest, str):
        if file_exists(dest):
            if show_error:
                log.error('dest exists, {}'.format(dest))
        elif os.path.exists(os.path.dirname(dest)):
            flag = True
        else:
            if show_error:
                log.error('dest not valid, {}'.format(dest))
    else:
        if show_error:
            log.error('dest expect str, got {}'.format(dest))
    # do-the-thing
    out = None
    if flag:
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        dest = os.path.abspath(os.path.expanduser(os.path.expandvars(dest)))
        # print('!A-1', dest)
        is_gzipped = file_is_gzipped(src)
        if decompress:
            if is_gzipped:
                with xopen(src, 'rb') as r, \
                    xopen(dest, 'wb') as w:
                    shutil.copyfileobj(r, w)
                out = dest
            else:
                if show_error:
                    log.error('src is not gzipped, skipped. {}'.format(src))
        else:
            if is_gzipped:
                if show_error:
                    log.error('src is gzipped, skipped. {}'.format(src))
            else:
                with xopen(src, 'rb') as r, \
                    xopen(dest, 'wb', threads=threads,
                          compresslevel=compresslevel) as w:
                    shutil.copyfileobj(r, w)
                out = dest
    return out


def read_file(x, nrows=0, skip=0, strip_white=True, comment=''):
    """
    Read plain text file
    save each line as list()
    Parameters
    ----------
    x:  str
        Path to a file
    nrows:  int
        The maximum number of rows to read
        default: [0], ignored
    skip:  int
        The number of lines of the data to skip before beginning to read data
        default: [0]
    strip_white:  bool
        Stripping of leading and trailing white space, default: [True]
    comment:  str
        A string of one character, default [''], empty
    """
    out = None
    if isinstance(x, str):
        if os.path.exists(x):
            if os.path.isdir(x):
                log.error('read_lines() failed, exptect <file>, got <dir>')
            else:
                l = []
                try:
                    i = 0
                    with open(x) as r:
                        for line in r:
                            i += 1
                            # skip rows
                            if skip > 0 and i <= skip:
                                continue
                            # white-spaces
                            s = line.strip()
                            if strip_white:
                                s = s.lstrip()
                            # comment
                            if len(comment) == 1 and s.startswith(comment):
                                continue
                            # nrows
                            if nrows > 0 and i > nrows:
                                break
                            # save to output
                            l.append(s)
                except IOError as e:
                    log.error(e)
                out = l
        else:
            log.error('read_lines() failed, file not exists')
    else:
        log.error('read_lins() failed, expect str, got {}'.format(
            type(x).__name__))
    return out


def fix_out_dir(x):
    """
    fix out_dir, if not "str", set "cwd()"
    """
    if not isinstance(x, str):
        x = pathlib.Path.cwd()
    x = os.path.abspath(x)
    if not os.path.exists(x):
        os.makedirs(x)
    return x


## 2. file info
def file_is_gzipped(x):
    """
    Check if the file is gzipped or not
    see answer: https://stackoverflow.com/a/3703300
    and on wiki: https://en.wikipedia.org/wiki/Gzip
    see also this blog: https://www.thinbug.com/q/3703276
    check the magic number for gzipped file: '1f8b'
    """
    if file_exists(x):
        with open(x, 'rb') as r:
            out = binascii.hexlify(r.read(2)) == b'1f8b'
    else:
        out = False
    return out


def file_prefix(x, with_dir=False):
    """
    Extract the prefix of file,
    compatible for None
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    remove extensions
    .gz, .fq.gz
    """
    if isinstance(x, str):
        if x.endswith('.gz') or x.endswith('.bz2'):
            x = os.path.splitext(x)[0]
        out = os.path.splitext(x)[0]
        if not with_dir:
            out = os.path.basename(out)
    elif isinstance(x, list):
        out = [file_prefix(i, with_dir) for i in x]
    elif x is None:
        out = None
    else:
        log.error('unknown x, str,list,None expected, got {}'.format(
            type(x).__name__))
        out = None
    return out


def file_abspath(x):
    """
    Return the absolute path of file
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    """
    if x is None or x == 'None': # in case toml format?!
        out = None
    elif isinstance(x, str):
        out = os.path.abspath(os.path.expanduser(x))
    elif isinstance(x, list):
        out = [file_abspath(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(x).__name__))
        out = x
    return out


def file_exists(x):
    """
    Check if file exists or not
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    """
    if x is None:
        out = False
    elif isinstance(x, str):
        out = os.path.exists(x) # file/dir/link
#         out = os.path.isfile(x)
    elif isinstance(x, list):
        out = [file_exists(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(file).__name__))
        out = False
    return out


def file_nrows(x):
    """
    Count the file rows
    count '\n'
    from @glglgl on stackoverflow, modified
    https://stackoverflow.com/a/9631635/2530783
    """
    def blocks(files, size = 1024 * 1024):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    if file_exists(x):
        with xopen(x, 'rt') as r:
            out = sum(bl.count('\n') for bl in blocks(r))
    else:
        log.error('x, file not exists, {}'.format(x))
        out = None
    return out


## 3. search files
def list_dir(x, full_names=True, recursive=False, include_dirs=False):
    """
    List all the files within the path
    see: list.dirs() in R
    see answers on :https://stackoverflow.com/a/3207973
    Parameters
    ----------
    x : str
        List files (dirs) in x
    full_names : bool
        Return the fullname of the files/dirs
    recursive : bool
        List files/dirs recursively
    include_dirs : bool
        Return the dirs
    """
    out = []
    if isinstance(x, str):
        if os.path.isdir(x):
            n = 0 # levels
            for (root, d, f) in os.walk(x):
                dirs = [os.path.join(root, i) for i in d] if full_names else d
                files = [os.path.join(root, i) for i in f] if full_names else f
                out += files
                if include_dirs:
                    out += dirs
                if not recursive:
                    break # first level
        else:
            log.error('list_dir() skipped, x not a directory: {}'.format(x))
    else:
        log.error('list_dir() skipped, x expect str, got {}'.format(
            type(x).__name__))
    return sorted(out)


def list_file(path='.', pattern='*', full_names=True, recursive=False,
    include_dirs=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()
    see base::list.files() function in R
    Parameters
    ----------
    x : str
        List files (dirs) in x
    pattern : str
        see pattern of fnmatch
        pattern:
        *       matches everything
        ?       matches any single character
        [seq]   matches any character in seq
        [!seq]  matches any char not in seq
        An initial period in FILENAME is not special.
        Both FILENAME and PATTERN are first case-normalized
        if the operating system requires it.
        If you don't want this, use fnmatchcase(FILENAME, PATTERN).
    full_names : bool
        Return the fullname of the files/dirs
    recursive : bool
        List files/dirs recursively
    include_dirs : bool
        Return the dirs
    example:
    list_file('./', '*.fq')
    """
    files = list_dir(path, full_names, recursive, include_dirs)
    files = [f for f in files if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(files)


## 4. valid files

def is_valid_file(x, fun):
    """
    Check if x is valid file: BAM|bigWig|BED

    Parameters
    ----------
    x : str
        Path to the file
    
    fun : function
        The function, is_valid_bam, is_valid_bigwig, is_valid_bed
    """
    if isinstance(x, str):
        out = fun(x)
    elif isinstance(x, list):
        out = all(list(map(fun, x)))
    else:
        out = False
    return out


def is_valid_bam(x):
    """
    Check if x is valid BAM file

    Parameters
    ----------
    x : str
        Path to the BAM file
    """
    out = False
    if isinstance(x, str):
        if os.path.exists(x):
            try:
                s = pysam.AlignmentFile(x)
                out = True
            except ValueError as err:
                out = False
    return out


def is_valid_bigwig(x):
    """
    Check if x is bigWig file
    retrieve by keywords
    bw files:
    """
    try:
        with pyBigWig.open(x) as bw:
            out = bw.isBigWig()
    except RuntimeError as err:
        out = False
        log.error(err)
    return out


def is_valid_bed(x):
    """
    Check if x is valid BED file
    file name: .bed, exists
    """
    out = False
    if isinstance(x, str):
        out1 = os.path.exists(x)
        x_ext = os.path.splitext(x)[1]
        out2 = x_ext.lower() in ['.bed', '.bed6', '.bed12', '.narrowpeak']
        out = out1 and out2
    return out


