#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import utils

sequenceFile = sys.argv[1]
clusterInfoFile = sys.argv[2]

sequences = {header : seq for header, seq in utils.parse_fasta(sequenceFile)}
cluster, centroids, failbob = utils.parse_clusterFile(clusterInfoFile)

unclusteredFile = '.'.join(sequenceFile.split(".")[:-1])+"_UNCLUSTERED.fasta"

with open(sequenceFile+"TEST", 'w') as outputStreamCluster:
  with open(unclusteredFile, 'w') as outputStream:
    for clusterID, accIDs in cluster.items():
      if len(accIDs) == 1:
        accID = accIDs[0]
        outputStream.write(f">{accID}\n{sequences[accID]}\n")
      else:
        for accID in accIDs:
          if accID in centroids:
            outputStreamCluster.write(f">{accID}\n{sequences[accID]}\n")