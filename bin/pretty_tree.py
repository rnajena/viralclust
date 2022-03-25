#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
from ete3 import Tree, TreeStyle


t = Tree(sys.argv[1], format=0)

ts = TreeStyle()
ts.mode = "c"
ts.arc_start = 0 # 0 degrees = 3 o'clock
ts.arc_span = 359
ts.optimal_scale_level = "full"
ts.show_leaf_name = False
#ts.scale = 5


t.show(tree_style = ts)
