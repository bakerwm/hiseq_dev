#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Demo for trackhub, to create CompositeTrack

Stage:

1. Create simple compositeTrack (for few tracks)
2. Create two compositeTracks
3. Add subgroups for tracks (for large number of tracks)
4. Automatic recognize subtrack, subgroups, ...
"""


import os
import sys
import re
import trackhub
from hiseq.utils.utils import log, update_obj, get_date
from hiseq.utils.file import file_abspath



class TrackFile(object):
    """
    Extract information from file

    Parameters
    ----------
    s : str
        The path to a track file,

    subgroups : dict
        The subgroups structure in dict format, including dimX, dimY, ...

    colorDim : str
        Assign colors to the tracks based on which dimension, dimX, dimY, ...
        default: [dimX]

    colorPal : int or str
        choose the group of colors, 1 to N, or the name of fish
        candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
        'Centropyge_loricula'], default: [1]

    label_rm_list : list
        A list of strings, remove from the 'short_label', to make sure the
        length of short label less than 17 characters. default: []
    """
    def __init__(self, s, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.s = s
        self.init_args()


    def init_args(self):
        default_args = {
            'subgroups': {},
            'colorDim': 'dimX',
            'colorPal': 1, # color palette, option: 1, 2, 3
            'label_rm_list': [],
        }
        self = update_obj(self, default_args, force=False)
        self.fname, self.fext = os.path.splitext(os.path.basename(self.s))
        self.fname = self.sanitize(self.fname)
        self.label = self.sanitize(self.fname, self.label_rm_list)
        self.ftype = self.filetype()
        # self.fcolor = self.pick_color(self.fname)
        self.subgroup = self.get_subgroup() # blank dict
        self.color = self.get_color() # '255,0,0'
        self.track = self.signal_track() if self.ftype in \
            ['bigWig', 'bedGraph', 'bigNarrowPeak'] else self.region_track()


    def filetype(self):
        """
        Determine the type of file
        Supported format list:
        - bigWig: bigWig, bigwig, bw
        - bigBed: bigBed, bigbed, bb
        - bedGraph: bedGraph, bedgraph, bg
        - bigNarrowPeak: bigNarrowPeak, np
        """
        ext = self.fext.lstrip('.').lower()
        # predefined ftypes
        ftypes = {
            'bw': 'bigWig',
            'bigwig': 'bigWig',
            'bb': 'bigBed',
            'bigbed': 'bigBed',
            'bg': 'bedGraph',
            'bdg': 'bedGraph',
        }
        ft = ftypes.get(ext, None)
#         if self.s.endswith('.narrowPeak.bb'):
#             ft = 'bigNarrowPeak'
        return ft
    
    
    def is_track_file(self):
        """
        trackhub supported files:
        - bigWig
        - bigBed
        - bedGraph
        - bigNarrowPeak
        """
        return any([self.is_bw(), self.is_bb(), self.is_bg(), self.is_bigNarrowPeak])


    def is_bw(self):
        return self.ftype == 'bigWig'


    def is_bb(self):
        return self.ftype == 'bigBed'


    def is_bg(self):
        return self.ftype == 'bedGraph'

    
    def is_bigNarrowPeak(self):
        return self.ftype == 'bigNarrowPeak'
        

    def sanitize(self, s, remove_list=None):
        """
        Sanitize a string, for shortLabel, letters, numbers allowed
        only [A-Za-z0-9]

        Convert spaces to unserscore, remove other characters

        Parameters
        -----------
        s : str
            String to sanitize

        remove_list : list or None
            if list, remove the items from string; if None, skipped
        """
        # remove non-characters:
        s1 = re.sub('[^A-Za-z0-9_ ]', '', s)
        # replace spaces to underscore
        s1 = re.sub('\s+', '_', s1)
        # remove specific strings
        if isinstance(remove_list, list) and len(remove_list) > 0:
            for i in remove_list:
                if not isinstance(i, str):
                    continue
                s1 = re.sub(i, '', s1, flags=re.IGNORECASE)
        return s1


    def get_subgroup(self):
        """
        Assign subgroups for file
        Searching 'mapping()' of each subgroup in filename
        return the 1st match
        Format for subgroups
        {
            'dimX': {
                'name': 'stage',
                'label': 'stage',
                'mapping': {
                    '2h': '2hrs',
                    '4h': '4hrs'
                }
            },
            'dimY': {
                'name': 'gene',
                'label': 'Gene',
                'mapping': {
                    'geneA': 'GeneA',
                    'geneB': 'GeneA'
                }
            }
        }
        # output format:
        # subgroup = {'stage': '2h', 'gene': 'geneA'}
        """
        subgroup = {}
        if isinstance(self.subgroups, dict) and len(self.subgroups) > 0:
            for d, v in self.subgroups.items():
                # iterate, dimensions
                # dimX, dimY, dimA, ...
                d_name = v.get('name', None)
                d_mapping = v.get('mapping', None)
                for m, n in d_mapping.items():
                    # iterate, groups
                    # 2h, 4h, ... (in dimX['mapping'])
                    if str(m) in self.fname:
                        subgroup[d_name] = m
                        break # stop iteration
                # check mapping
                if subgroup.get(d_name, None) is None:
                    log.error('subgroup_list not found in file: {}, {}'.format(
                        self.fname, list(d_mapping.keys())))
                    print('!A-5', d_name, self.s)
                    raise ValueError('subgroups failed')
        # format:
        # subgroup = {'stage': '2h', 'gene': 'geneA'}
        return subgroup


    def fish_colors(self, colorPal=1, n=10):
        """
        Pick colors
        Parameters
        ----------
        colorPal : int or str
            choose the group of colors, 1 to N, or the name of fish
            candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
            'Centropyge_loricula'], default: [1]

        n : int
            number of colors return, default: 10

        # colors from fishualize R package
        # https://nschiett.github.io/fishualize/articles/overview_colors.html

        Colors from R package:
        > fishualize::fish_pal(option = "Scarus_hoefleri")(10)
         [1] "#D2372CFF" "#E25719FF" "#F0780BFF" "#F89814FF"
         [5] "#EFB52BFF" "#C6CB43FF" "#7AD45CFF" "#00CE7DFF"
         [9] "#00BCABFF" "#0499EAFF"

        > fishualize::fish_pal(option = "Gramma_loreto")(10)
         [1] "#020122FF" "#1E2085FF" "#4029CBFF"
         [4] "#6628EEFF" "#901CEDFF" "#B804CAFF"
         [7] "#D61693FF" "#E6445DFF" "#EE7A30FF"
        [10] "#F0BF0BFF"

        > fishualize::fish_pal(option = "Centropyge_loricula")(10)
         [1] "#8F1D1EFF" "#B30029FF" "#DF002AFF"
         [4] "#FF7D1AFF" "#FFBD17FF" "#E7BE5AFF"
         [7] "#988591FF" "#0043A0FF" "#001A72FF"
        [10] "#000000FF"
        """
        # pre-defined colors
        c1 = ['#D2372C', '#E25719', '#F0780B', '#F89814', '#EFB52B',
            '#C6CB43', '#7AD45C', '#00CE7D', '#00BCAB', '#0499EA']
        c2 = ['#020122', '#1E2085', '#4029CB', '#6628EE', '#901CED',
            '#B804CA', '#D61693', '#E6445D', '#EE7A30', '#F0BF0B']
        c3 = ['#8F1D1E', '#B30029', '#DF002A', '#FF7D1A', '#FFBD17',
            '#E7BE5A', '#988591', '#0043A0', '#001A72', '#000000']
        # RGB (10 colors, repeat twice, 20 items)
        color_d = {
            'Scarus_hoefleri': list(map(self.hex2rgb, c1*2)),
            'Gramma_loreto': list(map(self.hex2rgb, c1*2)),
            'Centropyge_loricula': list(map(self.hex2rgb, c1*2))
        }
        # get the fish_list
        fish_list = list(color_d.keys())
        # determine the fish (colorBy)
        if isinstance(colorPal, int):
            colorPal = colorPal - 1 # to 0-indexed
            if not colorPal in range(len(fish_list)):
                colorPal = 0
            fish = fish_list[colorPal]
        elif isinstance(colorPal, str):
            fish = colorPal if colorPal in color_d else 'Scarus_hoefleri'
        else:
            fish = 'Scarus_hoefleri'
        # output colors
        return color_d.get(fish, c1)[:n]


    def hex2rgb(self, h):
        """
        Convert Hex to RGB code
        see: https://stackoverflow.com/a/29643643
        Parameters
        ----------
        h : str
            String hex color value
        >>> hex2rgb('#D2372C')
        """
        # extract the first 6 characters, exclude '#'
        if isinstance(h, str):
            if h.startswith('#') and len(h) >= 7:
                # h = h.lstrip('#')[:6]
                pass
            else:
                raise ValueError('hex color expected, got {}'.format(h))
        else:
            raise ValueError('hex color, str expected, got {}'.format(
                type(h).__name__))

        return ','.join(map(str,
            [int(h[1:3], 16),
            int(h[3:5], 16),
            int(h[5:7], 16)]
            ))


    def get_color(self, colorDim=None):
        """
        Get color, based on subgroups (dimX, dimY, dimA)
        Parameter
        ---------
        colorDim : str
            The dimension, to assign colors, dimX, dimY, dimA, ...
        subgroups:
        {
            'dimX': {
                'name': 'stage',
                'label': 'stage',
                'mapping': {
                    '2h': '2hrs',
                    '4h': '4hrs'
                }
            },
            'dimY': {
                'name': 'gene',
                'label': 'Gene',
                'maping': {
                    'geneA': 'GeneA',
                    'geneB': 'GeneA'
                }
            }
        }
        """
        if colorDim is None:
            colorDim = self.colorDim
        try:
            if not colorDim in self.subgroups:
                colorDim = 'dimX'
        except:
            print('!A-1, file: {}'.format(self.s))
            colorDim = 'dimX'
        # self.subgroups is the subgroups structure {dimX: , dimY: , ...}
        # extract the mapping from the colorDim (dimX)
        ttt = self.subgroups.get(colorDim, {})
        if ttt is None:
            print('!A-2', ttt)
            sys.exit()
        g_mapping = self.subgroups.get(colorDim, {}).get('mapping', {})
        # assign colors to each groups in mapping
        g_list = list(g_mapping.keys()) # 2h, 4h, ...
        c_list = self.fish_colors(colorPal=self.colorPal, n=len(g_list))
        # build the color dict for the mapping
        # {2h: 'rgb', 4h: 'rgb'}
        c_dict = {k:v for k, v in zip(g_list, c_list)}
        # self.subgroup, is groups for the file
        # format: stage=2h gene=geneA
        # extract the name of colorDim
        # dimX.name == 'stage'
        colorDimName = self.subgroups.get(colorDim, {}).get('name', None)
        # extract the group for colors
        # g = '2h'
        g = self.subgroup.get(colorDimName)
        # assign color
        return c_dict.get(g, '0,0,255')


    def signal_track(self):
        if self.ftype:
            track = trackhub.Track(
                name=self.fname+'_signal',
                short_label=self.label,
                source=self.s,
                visibility='full',
                viewLimits='0:10',
                maxHeightPixels='8:40:128',
                subgroups=self.subgroup,
                color=self.color,
                tracktype=self.ftype)
        else:
            track = None
        return track


    def region_track(self):
        if self.ftype:
            track = trackhub.Track(
                name=self.fname+'_region',
                short_label=self.label,
                source=self.s,
                visibility='dense',
                subgroups=self.subgroup,
                color=self.color,
                tracktype=self.ftype)
        else:
            track = None
        return track


    def to_track(self):
        """
        Choose {signal|region} Track, based on the filetype
        - signal Track, bigWig, bedGraph
        - region Track, bigBed
        """
        if self.ftype in ['bigWig', 'bedGraph', 'bigBed']:
            return self.track
        else:
            raise ValueError('bigWig, bedGraph, bigBed expected, got{}'.format(
                self.s))
