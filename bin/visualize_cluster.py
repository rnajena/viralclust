#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import itertools

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns

import utils


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)



sequenceFile = sys.argv[1]
clusterInfoFile = sys.argv[2]
cluster, centroids, failbob = utils.parse_clusterFile(clusterInfoFile)

allSequences = { header : seq for header, seq in utils.parse_fasta(sequenceFile) }
if len(sys.argv) == 4:
  goiSequences = { header : seq for header, seq in utils.parse_fasta(sys.argv[3]) }
  allSequences.update(goiSequences)


k = 7
nucleotides = set(["A","C","G","T"])
allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(nucleotides, repeat=k))}

profiles = {}

for header, sequence in allSequences.items():
  profile = [0]*len(allKmers)
  for kmer in iter([sequence[start : start + k] for start in range(len(sequence) - k)]):
      try:
        profile[allKmers[kmer]] += 1
      except KeyError:
        continue
  kmerSum = sum(profile)
  try:
    #profile = list(map(lambda x: x/kmerSum, profile))          
    profiles[header] = profile
  except ZeroDivisionError:
    print(header + " skipped, due to too many N's.")



for clusterID, header in cluster.items():
  subMatrixProfile = [ profiles[x] for x in header ]
  try:
    pairwise_distances = pd.DataFrame(data=squareform(pdist(subMatrixProfile)), index=header, columns=header)
  except ValueError:
    print(clusterID, header, len(subMatrixProfile))
    continue
  averages = np.mean(pairwise_distances, axis=0)
  pairwise_distances['average'] = averages
  sns.color_palette("mako", as_cmap=True)  
  ax = sns.heatmap(pairwise_distances) 
  

  for y_lab, x_lab in zip(ax.get_yticklabels(), ax.get_xticklabels()):
    y_text = y_lab.get_text()
    x_text = x_lab.get_text()
    if y_text in centroids and y_text not in goiSequences:
      y_lab.set_weight('bold')
      y_lab.set_color('purple')
      y_lab.set_size(20)
      x_lab.set_weight('bold')
      x_lab.set_color('purple')
      x_lab.set_size(20)

  plt.savefig(f"test_{clusterID}.png")
  plt.clf()
