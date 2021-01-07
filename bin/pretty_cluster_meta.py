#!/usr/bin/env python3
# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
from collections import Counter


def prettyCounterPrint(l, delim=', '):
  return(delim.join([f"{k}: {v}" for k,v in dict(Counter(l).most_common()).items()]))

inputFile = sys.argv[1]

cluster = ''
speciesList = []
locationList = []
dateList = []

holyTable = {}

with open(inputFile, 'r') as inputStream:
  for line in inputStream:
    line = line.strip()
    if line.startswith("## C"):
      if cluster:
        print()
        print(f'### {cluster} ###')
        print(prettyCounterPrint(speciesList))
        print(prettyCounterPrint(locationList))
        print(prettyCounterPrint(dateList))
      cluster = line.lstrip('#').strip()
      speciesList = []
      locationList = []
      dateList = []
      holyTable[cluster] = (speciesList, locationList, dateList)
      continue

    if line.startswith('###') or not line.strip():
      continue
  
    speciesList.append(line.split(',')[3])
    locationList.append(line.split(',')[4])
    dateList.append(line.split(',')[5])
