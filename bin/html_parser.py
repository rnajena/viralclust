#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


from collections import defaultdict

def parse_fasta(filePath):
  """
  """

  with open(filePath, 'r') as inputStream:
    header = ''
    seq = ''

    for line in inputStream:
      if line.startswith(">"):
        if header:
          yield (header, seq)

        header = line.rstrip("\n").replace(':','_').replace(' ','_').strip(">_")
        seq = ''
      else:
        seq += line.rstrip("\n").upper().replace('U','T')
    yield (header, seq)

def parse_clusterFile(clusterFile):
  """
  """

  centroids = []
  gois = []
  failbob = []
  acc2cluster = {}
  cluster = defaultdict(set)
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
        if line.endswith("GOI"):
          gois.append(accID)
        if clusterNumber == -1:
          acc2cluster[accID] = clusterNumber
          failbob.append(accID)
          continue
        if line.endswith('*'):# or line.endswith("GOI"):
          centroids.append(accID)

        acc2cluster[accID] = clusterNumber
        cluster[clusterNumber].add(accID)

  for clusterID, seqs in cluster.items():
    if len(seqs) == 1:
      failbob.extend(seqs)
  return(acc2cluster, centroids, failbob, cluster, gois)


def parse_metaInfo(filePath):
  """
  """

  holyTable = {}
  cluster = ''
  
  with open(filePath, 'r') as inputStream:
    for line in inputStream:
      line = line.strip()
      if line.startswith("## C"):
        cluster = int(line.lstrip('#').strip().split(' ')[1])
        continue

      if line.startswith('###') or not line.strip():
        continue
      content = line.split(',')[:6] + [cluster]
      holyTable[content[0]] = content[1:]
  return (holyTable)