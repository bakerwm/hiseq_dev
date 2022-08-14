#!/usr/bin/env python

"""
Generate bigWig files using deeptools
## Why this script?
In order to generate normalized forward/reverse bigWig files
dUTP library,
# (forward) -–filterRNAstrand=forward keeps minus-strand reads, -f 16
# (reverse) -–filterRNAstrand=reverse keeps plus-strand reads, -F 16
scaleFactor
normalizeUsing 
extendReads
smoothLength 
...
"""

import os
import argparse
# from bam2bw import Bam2bw
from hiseq.bam2bw.bam2bw import Bam2bw
from hiseq.utils.utils import log, Config


# Bam2bw template
def make_bam2bw_config(**kwargs):
    # for Bam2bw
    args_bw = {
        'bam': None,
        'out_dir': None,
        'out_prefix': 'metaplot',
        'scaleFactor': 1.0,
        'normalizeUsing': 'None',
        'binSize': 100,
        'numberOfProcessors': 4,
        'blackListFileName': None,
        'genome': None,
        'effectiveGenomeSize': None,
        'overwrite': False,
        'strand_specific': False, # for bigWig files, matrix, ...
        'extendReads': None,
        'centerReads': False,
    }
    args_bw.update(kwargs)
    return args_bw


def show_help(x):
    msg = '\n'.join([
        '-'*80,
        '# 1. Generate a template config file',
        '$ python make_bw.py -t -c {}'.format(x),
        '# 2. Modify the values in YAML' ,
        '# Attentation to the following fields:',
        '  - bam_list: ',
        '  - out_dir: ',
        '  - out_prefix: ',
        '  - normalizeUsing: ',
        '  - binSize: ',
        '  - blackListFileName: ',
        '  - effectiveGenomeSize: ',        
        '# 3. Run the command again',
        '$ python make_bw.py -c {}'.format(x),
        '-'*80,
    ])
    print(msg)


def get_args():
    example = '\n'.join([
        'Example:',
        '# 1. Generate template config file',
        '$ python make_bw.py -c a.yaml -t',
        '# 2: Run program',
        '# modify the config file `a.yaml` according to your data',
        '$ python make_bw.py -c a.yaml',
    ])
    parser = argparse.ArgumentParser(
        prog='make_bw', description='make_bw', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', default=None, required=True,
        help='configs in .yaml file')
    parser.add_argument('-b', '--bam', dest='bam', nargs='+', default=None,
        help='bam files, default: [None]')
    parser.add_argument('-o', dest='out_dir', required=False, default=None,
        help='directory to save bigWig file')
    parser.add_argument('-t', '--get-template', dest='get_template', action='store_true',
        help='Generate the template arguments')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite the plot files')
    return parser


def main():
    args = {
        'bam': None,
        'config': None,
        'get_template': False,
        'overwrite': False,
    }
    args.update(vars(get_args().parse_args()))
    # make sure config.yaml file
    if not isinstance(args['config'], str):
        raise ValueError('config, expect str, got {}'.format(
            type(args['config']).__name__
        ))
    if args['get_template']:
        if os.path.exists(args['config']) and not args['overwrite']:
            log.info('could write to config, file exists: {}'.format(
                args['config']
            ))
        else:
            # dump_yaml(make_bam2bw_config(), args['config'])
            Config().dump(make_bam2bw_config(), args['config'])
        show_help(args['config'])
    else:
        # d0 = load_yaml(args['config'])
        d0 = Config().load(args['config'])
        d1 = make_bam2bw_config()
        d1.update(d0) # from yaml
        # remove None values, from cmd
        if args['bam'] is None:
            args.pop('bam')
        if args['out_dir'] is None:
            args.pop('out_dir')
        d1.update(args) # from cmd
        Bam2bw(**d1).run()


if __name__ == '__main__':
    main()

# EOF
