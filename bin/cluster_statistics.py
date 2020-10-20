#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
Usage:
  cluster_statistics.py [options] <treeFile> <seqFile> <clusterFile>

Options:
  -h, --help        Prints this help and exits.
  --ncbi NCBI       Connects to the NCBI nuccore and taxonomy database to get taxonomic information based of the
                    header line of each fasta record. Note: Recommended, if a GenBank ID is part of the header.
  --toolName NAME   Used as prefix for NCBI taxonomy summary.

"""

import sys
import dendropy
import re
import numpy as np
import time

import utils

from docopt import docopt
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict

# def get_clade_by_rank(rank, taxIDs, target_record):
#   scientificNames = {}

#   for x in target_record:
#     if x['Rank'] == rank:
#       scientificNames[x['TaxId']] = x['ScientificName']
#     else:
#       for level in x['LineageEx']:
#         if level['Rank'] == rank:
#           scientificNames[x['TaxId']] = level['ScientificName']
#           break
#   accID2Name = {accessionID : scientificNames[taxID] if not taxID == 'XXXXX' else "unclassified" for accessionID, taxID in taxIDs.items()}
#   return accID2Name

def retrieve_taxonomy(prefix, accID2desc):

  avgClusterPerSpecies = defaultdict(set)
  avgClusterPerGenus = defaultdict(set)

  with open(f'{prefix}_taxonomy_info.txt', 'w') as outputStream:
    for clusterID, accessionIDs in cluster.items():
      outputStream.write(f"Cluster: {int(clusterID)+1}\n")
      for acc in accessionIDs:
        if acc in accID2desc:
          description = accID2desc[acc]
          if clusterID in realCluster:
            avgClusterPerSpecies[description.split(',')[0]].add(clusterID)
            avgClusterPerGenus[description.split(',')[1]].add(clusterID)
          outputStream.write(f"{acc},{description}\n")
        else:
          outputStream.write(f"{acc},--,--,--\n")
      outputStream.write(f"####################\n")
    return(avgClusterPerSpecies, avgClusterPerGenus)

#   speciesPerCluster = 0
#   genusPerCluster = 0
#   clusterPerSpecies = defaultdict(set)
#   clusterPerGenus = defaultdict(set)
# outputStream.write(f"Cluster: {int(clusterID)+1}\n")
# for acc, genusName in accID2genus.items():
#   outputStream.write(f"{acc},{genusName.replace(' ','_')},{accID2species[acc].replace(' ','_')},{accID2GBdescription[acc]}\n")
# for acc in invalidIDs:
#   outputStream.write(f"{acc}, --, --, --\n")
# outputStream.write(f"####################\n")



  # avgClusterPerSpecies = sum([len(x) for x in clusterPerSpecies.values()])/len(clusterPerSpecies)
  # avgClusterPerGenus = sum([len(x) for x in clusterPerGenus.values()])/len(clusterPerGenus)
  # return(avgClusterPerSpecies, avgClusterPerGenus)

#########################################################################

args = docopt(__doc__)


treeFile = args['<treeFile>']
seqFile = args['<seqFile>']
clusterFile = args['<clusterFile>']
NCBI = bool(args['--ncbi'])
PREFIX = args['--toolName']

allSequences = {header : seq for header,seq in utils.parse_fasta(seqFile)}
cluster, centroids, failbob = utils.parse_clusterFile(clusterFile)

# HDBScan specific. -1 describe sequences that are not clustered at all.
#try:
#  cluster.pop('-1')
#except KeyError:
#  pass



if not failbob:
  failbob = [cluster for idx,cluster in cluster.items() if len(cluster) == 1]
realCluster = {idx : cluster for idx,cluster in cluster.items() if len(cluster) != 1 and idx != '-1'}
  # print(failbob)

tree = dendropy.Tree.get(path=treeFile, schema='newick')
dm = tree.phylogenetic_distance_matrix()
allDistances = dm.distances()

allCluster = np.array([len(cl) for _,cl in realCluster.items()])

if NCBI:
  import pickle
  with open(args['--ncbi'], 'rb') as inputStream:
    accID2desc = pickle.load(inputStream)

  (clusterPerSpecies, clusterPerGenus) = retrieve_taxonomy(PREFIX, accID2desc)
  #avgClusterPerSpecies = sum([len(x) for x in clusterPerSpecies.values()])/len(clusterPerSpecies)
  avgClusterPerSpecies = np.mean([len(x) for x in clusterPerSpecies.values()])
  #avgClusterPerGenus = sum([len(x) for x in clusterPerGenus.values()])/len(clusterPerGenus)
  avgClusterPerGenus = np.mean([len(x) for x in clusterPerGenus.values()])

else:
  avgClusterPerSpecies = avgClusterPerGenus = '--'


print(f"{len(allSequences)},{len(realCluster)},{np.min(allCluster)},{np.max(allCluster)},{np.mean(allCluster):.2f},{np.median(allCluster)},{np.mean(allDistances):.3f},{len(failbob)},{avgClusterPerSpecies:.2f},{avgClusterPerGenus:.2f}")

