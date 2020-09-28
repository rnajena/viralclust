#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import re
import os

import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict

import utils

clusterFile = sys.argv[1]
basename = os.path.basename(clusterFile)
basename = os.path.splitext(basename)[0]
#if dirname:
ornamentFile = basename + "_ornaments.map"
cssFile = basename + "_cluster_css.map"
#else:
#  ornamentFile = "ornaments.map"
#  cssFile = "cluster_css.map"

# centroids = []
# failbob = []
# cluster = defaultdict(list)
# clusterNumber = -1

# with open(clusterFile, 'r') as inputStream:
#   for line in inputStream:
#     line = line.rstrip()
#     if line.startswith('>'):
#       clusterNumber = line.split(' ')[1]
#       continue
#     else:
#       accID = line.split('>')[1].split(' ')[0]
#       if clusterNumber == -1:
#         failbob.append(accID)
#         continue
#       if line.endswith('*'):
#         centroids.append(accID)
#       cluster[clusterNumber].append(accID)

cluster, centroids, failbob = utils.parse_clusterFile(clusterFile)

cmap = matplotlib.cm.get_cmap('Spectral')
norm = matplotlib.colors.Normalize(vmin=1,vmax=len(cluster))
colorMap = list(map(matplotlib.colors.rgb2hex, [cmap(norm(x))[:3] for x in range(len(cluster))]))


with open(cssFile, 'w') as outputStream:
  for idx, (number, sequences) in enumerate(cluster.items()):
    if number == '-1':
      continue
    outputStream.write(f'"stroke-width:1.5; stroke:{colorMap[idx]}" Clade {" ".join(sequences)}\n')

with open(ornamentFile, 'w') as outputStream:
  outputStream.write(f'"<circle style=\'fill:cyan;stroke:blue\' r=\'3\'/>" I {" ".join(centroids)}\n')
  if failbob:
    outputStream.write(f'"<circle style=\'fill:red;stroke:red\' r=\'3\'/>" I {" ".join(failbob)}\n')