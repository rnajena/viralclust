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
print(accessionsInCluster)

clusterInTarget = [cl for cl, items in targetContent.items() if any([x[0] in accessionsInCluster for x in items])]
print(clusterInTarget)