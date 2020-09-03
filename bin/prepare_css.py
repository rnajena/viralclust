#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import sys
import re

clusterFile = sys.argv[1]
centriods = []
cluster = {}

#genbankACCRegex = re.compile(r'[^A-Z]([A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{2}[0-9]{8}|NC_[0-9]{6})[^0-9]')

with open(clusterFile, 'r') as inputStream:
  for line in inputStream:
    if line.startswith('>'):
      clusterNumber = line.split(' ')[1]
      continue
    #else:
      #accessionID = list(set(re.findall(genbankACCRegex, line)))

      
      

      