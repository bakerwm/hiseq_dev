#!/usr/bin/env python3
# -*- encoding: utf8 -*- 

"""
list hiseq dir, file

hiseq list -n name -t _rx [dir1, ...]
"""

import os
import sys
import argparse
from hiseq.utils.utils import get_date, Config
from hiseq.utils.hiseq_utils import list_hiseq_dir, list_hiseq_file, is_hiseq_dir
from hiseq.utils.file import list_dir, file_abspath


def list_dirs(x, **kwargs):
    """
    list the following directories:
    current directory
    sub_folder level-1
    """
    args = {'recursive': False}
    args.update(kwargs)
    if isinstance(x, str):
        if os.path.isdir(x):
            out = [i for i in list_dir(x, include_dirs=True, recursive=args['recursive']) if os.path.isdir(i)]
            out.append(x)
    elif isinstance(x, list):
        out = [j for i in x for j in list_dirs(i, **kwargs)]
    else:
        out = []
    return out


def hiseq_list(x, **kwargs):
    """
    list hiseq files, dirs
    """
    args = {
        'name': None,
        'hiseq_type': 'auto',
        'recursive': False,
        'abspath': False,
        'out_json': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # 1. folders within x
    dirs = list_dirs(x, **kwargs) # update, root/subdir
    # 2. search hiseq dirs
    dx = []
    for d in dirs:
        if isinstance(args['name'], str):
            out = list_hiseq_file(d, args['name'], args['hiseq_type'])
            if out is None:
                out = []
        else:
            try:
                out = list_hiseq_dir(d, args['hiseq_type'])
            except:
                out = []
        if isinstance(out, str):
            out = [out]
        dx.extend(out)
    if args['abspath']:
        dx = file_abspath(dx)
    # only hiseq_dirs
    dx = [i for i in dx if is_hiseq_dir(i, args['hiseq_type'])]
    dx = sorted(list(set(dx)))
    # show message
    # args:
    a_in = [x] if isinstance(x, str) else x
    msg = '\n'.join([
        '='*80,
        f'{"Program":>14}: hiseq_list',
        f'{"Date":>14}: {get_date()}',
        f'{"hiseq_type":>14s}: {args["hiseq_type"]}',
        f'{"hiseq_name":>14s}: {args["name"]}',
        f'{"recursive":>14}: {args["recursive"]}',
        f'{"abspath":>14}: {args["abspath"]}',
        f'{"input_dir":>14}: ',
        '\n'.join([f'{"-":>14} {i}' for i in a_in]),
        f'{"output":>14}: {len(dx)} hits',
        '\n'.join([f'{"-":>14} {i}' for i in dx]),
        '='*80,
    ])
    print(msg, file=sys.stderr)
    # save to file
    dd = {
        'hiseq_list': dx,
        'hiseq_type': args['hiseq_type'],
        'hiseq_name': args['name'],
        }
    if isinstance(args['out_json'], str):
        if os.path.exists(args['out_json']) and not args['overwrite']:
            print(f'file exists: {args["out_json"]}')
        else:
            Config().dump(dd, args['out_json'])
    return dx


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. list the hiseq dirs',
        '$ hiseq hiseq_list path',
        ' ',
        '2. list the hiseq files',
        '$ hiseq hiseq_list -n smp_name path'        
    ])
    parser = argparse.ArgumentParser(
        prog='hiseq_list',
        description='list hiseq dir/file',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('x', nargs='+',
        help='list of directories') # positional arg
    parser.add_argument('-o', dest='out_json', required=True,
        help='save output in json file')
    parser.add_argument('-n', '--name', default=None, 
        help='name of the hiseq file, if not return the directory')
    parser.add_argument('-t', '--hiseq-type', dest='hiseq_type', default='auto',
        help='the hiseq type, default: [auto]')
    parser.add_argument('-r', '--recursive', dest='recursive', action='store_true',
        help='search the directories recursively')
    parser.add_argument('-a', '--abspath', action='store_true',
        help='return the absolute path')
    parser.add_argument('-O', '--overwrite', dest='overwrite',
        action='store_true', help='overwrite exists output json file')
    return parser


def main():
    args = vars(get_args().parse_args())
    out = hiseq_list(args.pop('x', None), **args)


if __name__ == '__main__':
    main()

#