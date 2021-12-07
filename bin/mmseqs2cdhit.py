#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
from collections import defaultdict

import utils

inputFile = sys.argv[1]
fastaFile = sys.argv[2]

if len(sys.argv) == 4:
  goiFile = sys.argv[3]
  goiHeader = [header for header, _ in utils.parse_fasta(goiFile)]

outputFile = f"{sys.argv[1]}.clstr"

clusterInfo = defaultdict(list)
fastaContent = {header : seq for header, seq in utils.parse_fasta(fastaFile)}
with open(inputFile, 'r') as inputStream:
  centroids = set()
  for line in inputStream:
    currentArray = line.split()
    centroid = currentArray[0]
    sequence = currentArray[1]
    centroids.add(centroid)
    clusterInfo[centroid].append(sequence)

with open(outputFile, 'w') as outputStream:
  for idx, centroid in enumerate(centroids):
    seqInCluster = clusterInfo[centroid]
    outputStream.write(f">Cluster {idx}\n")
    for seqIdx, seq in enumerate(seqInCluster):
      outputStream.write(f"{seqIdx}\t{len(fastaContent[seq])}nt, >{seq}")
      if seq == centroid:
        outputStream.write(' *')
      else:
        outputStream.write(" at +/13.37%")
      if len(sys.argv) == 4 and seq.strip('>') in goiHeader:
        outputStream.write("  GOI")
      outputStream.write("\n")