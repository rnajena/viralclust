#!/usr/bin/env python3
# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
from collections import Counter
import utils

def prettyCounterPrint(l, delim=', '):
  return( '_'.join(filter(None, delim.join([f"{k}: {v}" for k,v in dict(Counter(l).most_common()).items()]).split('_'))) )

inputFile = sys.argv[1]

holyTable = utils.parse_metaInfo(inputFile)

for cluster, elements in holyTable.items():
  speciesList = [x[3] for x in elements]
  locationList = [x[4] for x in elements]
  dateList = [x[5] for x in elements]
  print()
  print(f'### Cluster: {cluster} ###')
  print(prettyCounterPrint(speciesList))
  print(prettyCounterPrint(locationList))
  print(prettyCounterPrint(dateList))