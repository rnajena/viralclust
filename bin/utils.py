#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


import sys
import re
from collections import defaultdict

genbankACCRegex = re.compile(r'[^A-Z]*([A-Z]{2}[0-9]{8}|[A-Z]{2}[0-9]{6}|[A-Z][0-9]{5}|NC_[0-9]{6})[^0-9]?')

def reverseComplement(sequence):
  """

  """
  comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
  return(''.join([comp.get(x,'N')  for x in sequence.upper()[::-1]]))



def parse_fasta(filePath):
  """
  """

  with open(filePath, 'r') as inputStream:
    header = ''
    seq = ''

    for line in inputStream:
      if line.startswith(">"):
        if header:
          # seq = seq
          if seq.count('N') / len(seq) < 0.1:
            yield (header, seq)

        header = line.rstrip("\n ").replace(':','_').replace(' ','_').replace('|','_').lstrip(">").rstrip('_')
        accessionID = re.findall(genbankACCRegex, header)
        if accessionID:
          header = accessionID[0]
        else:
          header = header[:30]
        seq = ''
      else:
        seq += line.rstrip("\n").upper().replace('U','T')

    #seq = seq
    if seq.count('N') / len(seq) < 0.1:
      yield (header, seq)

def parse_clusterFile(clusterFile):
  centroids = []
  failbob = []
  cluster = defaultdict(list)
  clusterNumber = -1
  with open(clusterFile, 'r') as inputStream:
    for line in inputStream:
      line = line.rstrip()
      if line.startswith('>'):
        clusterNumber = int(line.split(' ')[1])
        #continue
      else:
        #print(line)
        accID = line.split('>')[1].split(' ')[0]
        if clusterNumber == -1:
          failbob.append(accID)
          continue
        if line.endswith('*') or line.endswith("GOI"):
          centroids.append(accID)
        cluster[clusterNumber].append(accID)

  return(cluster, centroids, failbob)

def parse_metaInfo(filePath):

  holyTable = {}
  cluster = ''
  content = []

  with open(filePath, 'r') as inputStream:
    for line in inputStream:
      line = line.strip()
      if line.startswith("## C"):
        if cluster:
          holyTable[cluster] = content
          content = []
        cluster = int(line.lstrip('#').strip().split(' ')[1])
        continue

      if line.startswith('###') or not line.strip():
        continue
      content.append(line.split(','))
    holyTable[cluster] = content
  return (holyTable)