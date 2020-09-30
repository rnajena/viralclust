#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

from collections import Counter
from multiprocessing import Pool
import multiprocessing as mp
import itertools
import math
import numpy as np
import os
import subprocess


import scipy
import random

import umap.umap_ as umap
import hdbscan

class Clusterer(object):
  """
  """

  id2header = {}
  d_profiles = {}
  header2id = {}
  dim = 0
  matrix = np.empty(shape=(dim,dim))

  genomeOfInterest = ''
  goiHeader = []
  goi2Cluster = {}

  def __init__(self, logger, sequenceFile, k, proc, outdir, subCluster=False, goi=""):
    """
    """

    self.subCluster = subCluster
    self.sequenceFile = sequenceFile
    self.outdir = outdir
    if not self.subCluster:
      self.reducedSequences = f"{self.outdir}/reduced.fasta"
    else:

      self.reducedSequences = f"{self.outdir}/{os.path.basename(sequenceFile)}"

    self.k = k

    nucleotides = set(["A","C","G","T"])
    self.allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(nucleotides, repeat=self.k))}
    self.proc = proc
    self.d_sequences = {}
    self.centroids = []
    self.allCluster = []
    self.clusterlabel = []
    self.probabilities = []

    if goi:
      Clusterer.genomeOfInterest = goi

  def __parse_fasta(self, filePath, goi=False):
    """
    """

    with open(filePath, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
            yield (header, seq)

          header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
          seq = ''
        else:
          seq += line.rstrip("\n").upper().replace('U','T')
      yield (header, seq)

  def read_sequences(self):
    """
    """
    idHead = -1
    fastaContent = {}
    for header, sequence in self.__parse_fasta(self.sequenceFile):
      if not self.subCluster:
        idHead += 1
        Clusterer.id2header[idHead] = header
        Clusterer.header2id[header] = idHead
        fastaContent[idHead] = sequence
      else:
        fastaContent[Clusterer.header2id[header]] = sequence

    if Clusterer.genomeOfInterest and not self.subCluster:
      for header, sequence in self.__parse_fasta(Clusterer.genomeOfInterest):
          idHead += 1
          Clusterer.goiHeader.append(idHead)
          Clusterer.id2header[idHead] = header
          Clusterer.header2id[header] = idHead
          fastaContent[idHead] = sequence

    if not self.subCluster:
      Clusterer.dim = len(fastaContent)
      Clusterer.matrix = np.zeros(shape=(Clusterer.dim, Clusterer.dim), dtype=float)

    if len(fastaContent) < 21:
      return 1
    return fastaContent


  def profile(self, entry):
    """
    """
    header, sequence = entry
    profile = [0]*len(self.allKmers)
    for k in iter([sequence[start : start + self.k] for start in range(len(sequence) - self.k)]):
        try:
          profile[self.allKmers[k]] += 1
        except KeyError:
          continue
    kmerSum = sum(profile)
    profile = list(map(lambda x: x/kmerSum, profile))
    return (header, profile)

  def determine_profile(self, proc):

    self.d_sequences = self.read_sequences()
    if self.d_sequences == 1:
      import shutil
      shutil.copyfile(self.sequenceFile, f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta')
      return 1

    p = Pool(self.proc)
    allProfiles = p.map(self.profile, self.d_sequences.items())
    p.close()
    p.join()
    for header, profile in allProfiles:
      Clusterer.d_profiles[header] = profile

  def calc_pd(self, seqs):
    """
    """
    for element in seqs:
      try:
        stuff = (element[0], element[1])
      except TypeError:
        return None
    seq1, profile1 = seqs[0]
    seq2, profile2 = seqs[1]
    distance = scipy.spatial.distance.cosine(profile1, profile2)
    return (seq1, seq2, distance)
    #

  def apply_umap(self):
    """
    """
    profiles = [(idx,profile) for idx, profile in Clusterer.d_profiles.items() if idx in self.d_sequences]
    vector = [x[1] for x in profiles]

    if self.subCluster:
      neighbors, dist = 5, 0.0
    else:
      neighbors, dist = 50, 0.25


    try:
      clusterable_embedding = umap.UMAP(
            n_neighbors=neighbors,
            min_dist=dist,
            n_components=20,
            random_state=42,
            metric='cosine',
        ).fit_transform(vector)

      clusterer = hdbscan.HDBSCAN()
      clusterer.fit(clusterable_embedding)

      self.clusterlabel = clusterer.labels_
      self.probabilities = clusterer.probabilities_
      if len(set(self.clusterlabel)) == 1:
        raise TypeError

    except TypeError:
      import shutil
      shutil.copyfile(self.sequenceFile, f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta')
      return 1

    self.allCluster = list(zip([x[0] for x in profiles], self.clusterlabel))
    if not self.subCluster:
      with open(f'{self.outdir}/cluster.txt', 'w') as outStream:
        for i in set(self.clusterlabel):
          with open(f'{self.outdir}/cluster{i}.fasta', 'w') as fastaOut:
            outStream.write(f">Cluster {i}\n")
            for idx, label in self.allCluster:
              if label == i:
                if idx in Clusterer.goiHeader:
                  Clusterer.goi2Cluster[Clusterer.id2header[idx]] = i
                outStream.write(f"{Clusterer.id2header[idx]}\n")
                fastaOut.write(f">{Clusterer.id2header[idx]}\n{self.d_sequences[idx].split('X'*10)[0]}\n")
          outStream.write("\n")


  def get_centroids(self, proc):
    """
    """
    seqCluster = { x : [] for x in set(self.clusterlabel)}

    for idx, cluster in self.allCluster:
      seqCluster[cluster].append(idx)

    p = Pool(self.proc)

    for cluster, sequences in seqCluster.items():
      if cluster == -1:
        for sequence in sequences:
          if sequence in Clusterer.goiHeader:
            self.centroids.append(sequence)
        continue #print(sequence)

      subProfiles = {seq : profile for seq, profile in Clusterer.d_profiles.items() if seq in sequences}

      if not self.subCluster:
        for result in p.map(self.calc_pd, itertools.combinations(subProfiles.items(), 2)):
          seq1, seq2, dist = result
          Clusterer.matrix[seq1][seq2] = dist
          Clusterer.matrix[seq2][seq1] = dist

      tmpMinimum = math.inf
      centroidOfCluster = -1

      if len(sequences) == 1:
        centroidOfCluster = cluster[0]
        self.centroids.append(centroidOfCluster)
        continue

      for sequence in sequences:
        averagedDistance = 0

        if sequence in Clusterer.goiHeader:
          #print(Clusterer.goiHeader)
          self.centroids.append(sequence)
          continue

        for neighborSequence in sequences:
          if sequence == neighborSequence:
            continue
          averagedDistance += Clusterer.matrix[sequence][neighborSequence]


        averagedDistance /= len(sequences)-1

        if averagedDistance < tmpMinimum:
          tmpMinimum = averagedDistance
          centroidOfCluster = sequence

      self.centroids.append(centroidOfCluster)
    p.close()
    p.join()


  def output_centroids(self):
    """
    """
    reprSeqs = {centroid : self.d_sequences[centroid] for centroid in self.centroids}
    outputPath = f'{self.outdir}/{os.path.splitext(os.path.basename(self.sequenceFile))[0]}_hdbscan.fasta'

    with open(outputPath, 'w') as outStream:
      for centroidID, sequence in reprSeqs.items():
        outStream.write(f">{Clusterer.id2header[centroidID]}\n{sequence}\n")

    with open(f'{outputPath}.clstr', 'w') as outStream:
      for i in set(self.clusterlabel):
        outStream.write(f">Cluster {i}\n")
        seqInCluster = [idx for idx,label in self.allCluster if label == i]
        for idx, seqIdx in enumerate(seqInCluster):
          outStream.write(f"{idx}\t{len(self.d_sequences[seqIdx])}nt, >{Clusterer.id2header[seqIdx]} ")
          if seqIdx in self.centroids:
            outStream.write('*\n')
          else:
            outStream.write("at +/13.37%\n")