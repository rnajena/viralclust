#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import utils

clusterSource = sys.argv[1]
clusterNumber = int(sys.argv[2])
clusterTarget = sys.argv[3]

sourceContent = utils.parse_metaInfo(clusterSource)
targetContent = utils.parse_metaInfo(clusterTarget)

accessionsInCluster = [x[0] for x in sourceContent[clusterNumber]]
clusterInTarget = [cl for cl, items in targetContent.items() if any([x[0] in accessionsInCluster for x in items])]
targetGenomes = [y[0] for x in clusterInTarget for y in targetContent[x]]
accessionsInTarget = len(targetGenomes)
bidirectionalAcc = [cl for cl, items in sourceContent.items() if any([x[0] in targetGenomes for x in items])]

print()
print(f"Genomes clustered in cluster {clusterNumber}: {len(accessionsInCluster)}\n")
print(f"These genomes are found in {clusterTarget} in the following cluster: {', '.join(map(str, clusterInTarget))}, summing up to {accessionsInTarget} genomes")
print(f"In turn, these {accessionsInTarget} genomes are found in cluster {', '.join(map(str, bidirectionalAcc))} in {clusterSource}")
