#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import dendropy

import utils

treeFile = sys.argv[1]
seqFile = sys.argv[2]
centroidFile = sys.argv[3]
clusterFile = sys.argv[4]


allSequences = {header : seq for header,seq in utils.parse_fasta(seqFile)}
centroidSequences = {header : seq for header,seq in utils.parse_fasta(centroidFile)}
cluster, centroids, failbob = utils.parse_clusterFile(clusterFile)

# HDBScan specific. -1 describe sequences that are not clustered at all.
try:
  cluster.pop('-1')
except KeyError:
  pass


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
    
print(len(allSequences), len(centroids), avgOverallSum)