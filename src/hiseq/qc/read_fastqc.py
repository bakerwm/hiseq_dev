#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Parse the plain text output of 
1. Falco [https://github.com/smithlabcode/falco]
2. FastQC []
file: fastqc_data.txt 
convert to Json format
Date: 2022-08-31
"""


# import os
# import sys
# import re
import json

# import plotly.express as px
# import pandas as pd
# from hiseq.utils.utils import Config, log, run_shell_cmd


class ReadFastQC(object):
    def __init__(self, x):
        self.x = x
        self.is_falco = self.is_valid_falco()
        self.is_fastqc = self.is_valid_fastqc()
        self.is_valid = self.is_falco or self.is_fastqc
        self.d = self.parse_fastqc(x) if self.is_valid else None

    def is_valid_falco(self):
        try:
            with open(self.x) as r:
                s = r.readline()  # 1st line
            out = s.startswith("##Falco")
        except:
            out = False
        return out

    def is_valid_fastqc(self):
        try:
            with open(self.x) as r:
                s = r.readline()  # 1st line
            out = s.startswith("##FastQC")
        except:
            out = False
        return out

    def tab_to_dict2(self, x):
        """
        Parameters:
        -----------
        x : list
            list of table, each row as a item in list
        Convert tab to dict
        column-1: key
        column-2: value
        """
        if not isinstance(x, list):
            return None
        if len(x) < 1:
            return None
        # check header
        if x[0].startswith("#"):
            x = x[1:]
        # split values
        d = {}
        for s in x:
            k, v = s.strip().split("\t", 1)
            d.update({k: v})
        return d

    def tab_to_dict(self, x):
        """
        Parameters:
        -----------
        x : list
            list of table, each row as a item in list
        Convert tab to dict
        #header1 header2 ...
        val1 val2 ...
        """
        if not isinstance(x, list):
            return None
        if len(x) < 1:
            return None
        # check header
        if x[0].startswith("#"):
            x1 = x[0][1:].strip()
            h = x1.split("\t")
            v = x[1:]
        else:
            # X0, X1, ...
            tabs = x[0].split("\t")
            h = ["X{}".format(i + 1) for i, _ in enumerate(tabs)]
            v = x
        # split values
        d = {}
        for s in v:
            tabs = s.strip().split("\t")
            # check, if h-size == tabs-size
            for i, j in enumerate(tabs):
                try:
                    j2 = d.get(h[i], [])
                except:
                    # print('!A-1', j2)
                    pass
                j2.append(j)
                d.update({h[i]: j2})
        return d

    def parse_fastqc(self, x):
        """
        Parameters:
        -----------
        x : str
            file of Falco output, eg: fastqc_data.txt
        The text output of Falco, FastQC
        Handle exceptions?
        """
        d = {}  # final dict
        with open(x) as r:
            v = []  # save table
            d2 = None  # for "Total Deduplicated Percentage" only
            for l in r:
                l = l.strip()
                if l.startswith("##"):
                    # 1st line: ##Falco	0.2.4
                    tool, ver = l[2:].split("\t")
                    d.update({"software": tool, "version": ver})
                    continue  # skip
                if l.startswith(">>"):
                    if l.startswith(">>END_MODULE"):
                        t2d = (
                            self.tab_to_dict2
                            if m == "Basic Statistics"
                            else self.tab_to_dict
                        )
                        d1 = {"module": m, "status": s, "table": t2d(v)}
                        if isinstance(d2, dict):
                            d1.update(d2)
                        d.update({m: d1})
                        v = []  # init d2
                        d2 = None  # init d2
                    else:
                        m, s = l[2:].split("\t", 1)  # update module/status
                    continue
                # for "Sequence Duplication Levels"
                if l.startswith("#Total Deduplicated Percentage"):
                    a, b = l[1:].strip().split("\t")
                    d2 = {a: b}
                else:
                    v.append(l)
        return d

    def list_modules(self):
        if isinstance(self.d, dict):
            [self.d.pop(i, None) for i in ["software", "version"]]
            return list(self.d.keys())

    def get_module(self, m):
        if isinstance(self.d, dict):
            return self.d.get(m, None)

    def to_json(self):
        if isinstance(self.d, dict):
            return json.dumps(self.d, indent=4)

    def save_as(self, out):
        try:
            with open(out, "wt") as w:
                json.dump(self.d, w, indent=4, sort_keys=False)
        except:
            print("Could not write to file: {}".format(out))


def main():
    x = "PROseq_01A_rep1_fastqc/fastqc_data.txt"
    m = "Per base sequence content"
    aa = ReadFastQC(x)
    print(aa.to_json())


if __name__ == "__main__":
    main()

#
