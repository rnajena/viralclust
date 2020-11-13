#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import utils

inputFile = sys.argv[1]
outputFile = f"{sys.argv[1]}.clstr"

if len(sys.argv) == 3:
  goiFile = sys.argv[2]
  goiHeader = [header for header, _ in utils.parse_fasta(goiFile)]

centroids = {}
cluster = {}

with open(inputFile, 'r') as inputStream:
  for line in inputStream:
    lineArray = line.strip().split('\t')
    if line.startswith('S'):
      centroids[lineArray[1]] = (lineArray[8], lineArray[2])

    if line.startswith('H'):
      if lineArray[1] in cluster:
        cluster[lineArray[1]].append((lineArray[8], lineArray[2], lineArray[3]))
      else:
        cluster[lineArray[1]] = [(lineArray[8], lineArray[2], lineArray[3])]

with open(outputFile, 'w') as outputStream:
  for idx, content in centroids.items():
    if len(sys.argv) == 3 and content[0] in goiHeader:
      goiString = "  GOI"
    else:
      goiString = ''
    outputStream.write(f">Cluster {idx}\n")
    outputStream.write(f"0\t{content[1]}nt, >{content[0]} *{goiString}\n")
    if not idx in cluster:
      continue

    for seqIdx, clustercontent in enumerate(cluster[idx]):
      if len(sys.argv) == 3 and clustercontent[0] in goiHeader:
        goiString = "  GOI"
      else:
        goiString = ''
      outputStream.write(f"{seqIdx+1}\t{clustercontent[1]}nt, >{clustercontent[0]} at +/{clustercontent[2]}%{goiString}\n")
