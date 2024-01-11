#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
CnR pipeline: level-1 (main port)
mission-1: generate design.toml
mission-2: run_pipe, parsing config from design.toml
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.utils.hiseq_utils import HiSeqDesignCnr
from hiseq.utils.file import fix_out_dir
from hiseq.cnr.cnr_rx import CnrRx
from hiseq.cnr.cnr_args import get_args_cnr
from hiseq.utils.utils import Config, update_obj, log


class Cnr(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)

    def run_single_rx(self, x):
        g = self.fq_groups.get(x, None)  # input, ip
        if not isinstance(g, dict):
            log.error("error in design: {}".format(x))
            return None
        # update ip, input
        args = self.__dict__.copy()
        # ip
        g_ip = g.get("ip")
        if isinstance(g_ip, dict):
            args.update(
                {
                    "build_design": False,
                    "ip_fq1": [v[0] for k, v in g_ip.items()],
                    "ip_fq2": [v[1] for k, v in g_ip.items()],
                }
            )
        g_input = g.get("input")
        if isinstance(g_input, dict):
            args.update(
                {
                    "input_fq1": [v[0] for k, v in g_input.items()],
                    "input_fq2": [v[1] for k, v in g_input.items()],
                }
            )
        if len(self.fq_groups) > 1:
            args["parallel_jobs"] = 1  # force
        CnrRx(**args).run()

    def run(self):
        if self.build_design:
            HiSeqDesignCnr(**self.__dict__).run()
        else:
            self.fq_groups = Config().load(self.design)
            if self.parallel_jobs > 1 and len(self.fq_groups) > 1:
                with Pool(processes=self.parallel_jobs) as pool:
                    pool.map(self.run_single_rx, self.fq_groups)
            elif len(self.fq_groups) >= 1:
                for i in list(self.fq_groups):
                    self.run_single_rx(i)
            else:
                raise ValueError("no data found: {}".format(self.design))


class CnrConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_init = {
            "build_design": False,
            "design": None,  # required
            "fq_dir": None,
            "ip": None,  # str
            "input": None,  # str
            # 'threads': 1,
            "parallel_jobs": 1,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = "cnr_ra"
        self.out_dir = fix_out_dir(self.out_dir)
        if not isinstance(self.design, str):
            raise ValueError(
                "--design, expect str, got {}".format(
                    type(self.design).__name__
                )
            )


def get_args():
    return get_args_cnr()


def main():
    args = vars(get_args().parse_args())
    Cnr(**args).run()


if __name__ == "__main__":
    main()

#
