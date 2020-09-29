#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import dendropy
import numpy as np

import utils

treeFile = sys.argv[1]
seqFile = sys.argv[2]
#centroidFile = sys.argv[2]
clusterFile = sys.argv[3]


allSequences = {header : seq for header,seq in utils.parse_fasta(seqFile)}
#centroidSequences = {header : seq for header,seq in utils.parse_fasta(centroidFile)}
cluster, centroids, failbob = utils.parse_clusterFile(clusterFile)

# HDBScan specific. -1 describe sequences that are not clustered at all.
try:
  cluster.pop('-1')
except KeyError:
  pass

if not failbob:
  failbob = [cluster[0] for _,cluster in cluster.items() if len(cluster) == 1]
  cluster = {idx : cluster for idx,cluster in cluster.items() if len(cluster) != 1}
  # print(failbob)

tree = dendropy.Tree.get(path=treeFile, schema='newick')
dm = tree.phylogenetic_distance_matrix()

#print(tree.leaf_nodes())

overallSum = 0
for idx, cl in cluster.items():
  clusterSum = 0
  for header in cl:
    if header in centroids:
      centroid = header
      centroid_taxon = tree.find_node_with_taxon_label(centroid.replace('_',' ')).taxon
      break
  for header in cl:
    if header == centroid:
      continue
    t1 = tree.find_node_with_taxon_label(header.replace('_',' ')).taxon
    clusterSum += (dm.distance(t1,centroid_taxon))
  overallSum += clusterSum / len(cl)
avgOverallSum = overallSum / len(cluster)

allCluster = np.array([len(cl) for _,cl in cluster.items()])
# print(allCluster)

print(f"{len(allSequences)}, {len(cluster)}, {np.min(allCluster)}, {np.max(allCluster)}, {np.mean(allCluster)}, {np.median(allCluster)}, {avgOverallSum}, {len(failbob)}")