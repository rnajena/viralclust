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

fastaContent = {}
clusterInfo = {}

with open(inputFile, 'r') as inputStream:
  header = ''
  seq = ''
  for line in inputStream:
    if line.startswith('>'):
      if not header:
        header = line.rstrip()
      else:
        fastaContent[header.split(' ')[0]] = seq
        seq = ''
        header = line.rstrip()
      headerArray = header.split(';')

      clusterInfo[headerArray[0].split(' ')[0]] = {y[0].strip() : y[1].strip() for y in [x.split('=') for x in headerArray[1:] if x]}
    else:
      seq += line.rstrip()
  fastaContent[header.split(' ')[0]] = seq

centroids = [header for header, headerInfo in clusterInfo.items() if headerInfo['cluster_center'] == 'True']

with open(outputFile, 'w') as outputStream:
  for idx, centroid in enumerate(centroids):
    seqInCluster = [header for header, headerInfo in clusterInfo.items() if headerInfo['cluster'] == centroid.lstrip('>')]
    outputStream.write(f">Cluster {idx}\n")
    for seqIdx, header in enumerate(seqInCluster):
      outputStream.write(f"{seqIdx}\t{len(fastaContent[header])}nt, {header}")
      if header == centroid:
        outputStream.write(' *')
      else:
        outputStream.write(" at +/13.37%")
      if len(sys.argv) == 3 and header.strip('>') in goiHeader:
        outputStream.write("  GOI")
      outputStream.write("\n")