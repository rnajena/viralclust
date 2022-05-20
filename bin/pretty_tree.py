#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
from ete3 import PhyloTree, TreeStyle, NodeStyle


t = PhyloTree(sys.argv[1], format=0)

ts = TreeStyle()
ts.mode = "c"
ts.arc_start = 180 # 0 degrees = 3 o'clock
ts.arc_span = 359
ts.optimal_scale_level = "full"
ts.root_opening_factor = 0.01
ts.show_leaf_name = False

nstyle = NodeStyle()
nstyle["size"] = 0

for node in t.traverse():
  node.set_style(nstyle)

t.unroot()
t.show(tree_style = ts)
