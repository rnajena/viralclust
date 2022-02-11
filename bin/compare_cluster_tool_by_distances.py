#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import itertools
import math
import tempfile

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator


import utils

sequenceFile = sys.argv[1]
clusterInfoFile = sys.argv[2]
cluster, centroids, failbob = utils.parse_clusterFile(clusterInfoFile)

allSequences = { header : seq for header, seq in utils.parse_fasta(sequenceFile) }
if len(sys.argv) == 4:
  goiSequences = { header : seq for header, seq in utils.parse_fasta(sys.argv[3]) }
  allSequences.update(goiSequences)

for clusterID, header in cluster.items():
  if(len(header) == 1):
    continue
  
  pseudoFasta = ''.join([f">{x}\n{seq}\n" for x,seq in allSequences.items() if x in header])
  tmpFile = tempfile.NamedTemporaryFile()
  with open(tmpFile.name, 'w') as outputStream:
    outputStream.write(pseudoFasta)

  #print("Doing the alignment")
  mafft_cline = MafftCommandline(input=tmpFile.name, treeout=False, thread=4)
  stdout, stderr = mafft_cline()
  tmpAlignment = tempfile.NamedTemporaryFile()
  with open(tmpAlignment.name, 'w') as outputStream:
    outputStream.write(stdout)

  align = AlignIO.read(tmpAlignment.name, "fasta")
  calc = DistanceCalculator()
  #print("Calculating the distances")
  dm = calc.get_distance(align)

  pairwise_distances = pd.DataFrame(dm.matrix, index=dm.names, columns=dm.names).fillna(0)
  pairwise_distances = pairwise_distances + pairwise_distances.T

  averages = np.mean(pairwise_distances, axis=0)
  pairwise_distances['average'] = averages
  minDistanceSeq = pairwise_distances['average'].idxmin()

  for i,row in pairwise_distances.iterrows():
    #print(i, row['average'], end=' ')
    if i == minDistanceSeq:
      minimum = row['average']
    if i in centroids:
      centroid = row['average']
      centroidSeq = i
    if i in goiSequences:
      goi = row['average']
      goiSeq = i

  print(f"{clusterID},{len(header)},{np.mean(row['average'])},{minimum},{centroid},{minDistanceSeq==centroidSeq},{np.absolute(minimum-centroid)}")
