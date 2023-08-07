#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
ViralClust is a python program that takes several genome sequences
from different viruses as an input.
These will be clustered these sequences into groups (clades) based
on their sequence similarity. For each clade, the centroid sequence is
determined as representative genome, i.e. the sequence with the lowest
distance to all other sequences of this clade.

Python Dependencies:
  docopt
  BioPython
  colorlog
  numpy
  scipy
  umap-learn
  hdbscan

Contact:
  kevin.lamkiewicz@uni-jena.de


Usage:
  hdbscan_virus.py [options] <inputSequences> [<genomeOfInterest>]

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
  --version                               Prints the version of viralClust and exits.
  -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]
  --metric METRIC                         Distance metric applied by UMAP (if applied) and HDBSCAN.
                                          The following are supported:
                                          'euclidean', 'manhattan', 'chebyshev', 'minkwoski',
                                          'canberra', 'braycurtis',
                                          'mahalanobis', 'seuclidean',
                                          'cosine'.
                                          If an invalid metric is set, ViralClust will default back to 
                                          the cosine distance.
                                          [Default: cosine]

  --umap                                  Flag that determines whether (instead of PCA) a UMAP analysis is used for dimension reduction. [Default: False]

  --neighbors NEIGHBORS                   Number of neighbors considered by UMAP to reduce the dimension space.
                                          Low numbers here mean focus on local structures within the data, whereas 
                                          larger numbers may loose fine details. [default: 50]
  --dThreshold dTHRESHOLD                 Sets the threshold for the minimum distance of two points in the low-dimensional space.
                                          Smaller thresholds are recommended for clustering and identifying finer topological structures
                                          in the data. [Default: 0.25]
  --dimension DIMENSION                   UMAP tries to find an embedding for the input data that can be represented by a low-dimensional space.
                                          This parameter tells UMAP how many dimensions should be used for the embedding. Lower numbers may result 
                                          in loss of information, whereas larger numbers will increase the runtime. [Default: 20]

  --clusterSize CLUSTERSIZE               This parameter forces HDBSCAN to form cluster with a size larger-equal to CLUSTERSIZE.
                                          Be aware that some data points (i.e. genomes) or even whole subcluster are considered as noise, if this parameter is set too high.
                                          E.g., if a very distinct viral genus has 40 genomes and the parameter is set to anything >40, HDBSCAN will not form
                                          the genus specific cluster. [Default: 5]
  --minSample MINSAMPLE                   Intuitively, this parameter declares how conservative clustering is performed. Higher values will lead 
                                          to more points considered noise, whereas a low value causes "edge-cases" to be grouped into a cluster.
                                          The default parameter is the same as CLUSTERSIZE. [Default: CLUSTERSIZE]
  
  
  

  

"""

import sys
import os
import logging
import glob
import shutil
import csv
import time
import gc
import pickle

import numpy as np
from multiprocessing import Pool

from datetime import datetime

from colorlog import ColoredFormatter
from docopt import docopt
from Bio import Phylo

from collections import Counter
import multiprocessing as mp
import itertools
import math
import subprocess

import scipy
import random

import umap.umap_ as umap
import hdbscan
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA

def main():
  """
  """

  global logger
  global scipyDistances
  global id2header
  global header2id
  global goiHeader
  global goi2Cluster
  global metric

  logger = __create_logger()
  (inputSequences, goi, outdir,  KMER, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, umap_flag) = __parse_arguments(docopt(__doc__))

  id2header = {}
  header2id = {}
  goiHeader = []
  goi2Cluster = {}

  scipyDistances = {
    'euclidean' : scipy.spatial.distance.euclidean ,
    'manhattan' : scipy.spatial.distance.cityblock ,
    'chebyshev' : scipy.spatial.distance.chebyshev  ,
    'minkwoski': scipy.spatial.distance.minkowski ,
    'canberra' : scipy.spatial.distance.canberra ,
    'braycurtis' : scipy.spatial.distance.braycurtis ,
    'mahalanobis' : scipy.spatial.distance.mahalanobis ,
    'seuclidean' : scipy.spatial.distance.seuclidean ,
    'cosine' : scipy.spatial.distance.cosine  
  }

  reducedSequences = f"{outdir}/{os.path.basename(inputSequences)}"
  nucleotides = set(["A","C","G","T"])
  allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(nucleotides, repeat=KMER))}

  logger.info("Starting to cluster you data. Stay tuned.")

  logger.info("Determining k-mer profiles for all sequences.")
  numberOfSequences = sequences_to_profiles(inputSequences, outdir, allKmers, KMER)

  if goi:
    logger.info("Determining k-mer profiles for genomes of interest.")
    numberOfSequences = sequences_to_profiles(goi, outdir, allKmers, KMER, previousID=numberOfSequences)
  
  if umap_flag:
    logger.info("Clustering with UMAP and HDBSCAN.")
  else:
    logger.info("Clustering with PCA and HDBSCAN.")
  embedding = dimension_reduction(outdir, neighbors, threshold, metric, umap_flag)
  clusteredSeqs = cluster_with_hdbscan(embedding, metric)

  logger.info(f"The {len(id2header)} sequences were clustered into {len(set([x[1] for x in clusteredSeqs]))} groups.")
  logger.info(f"There are {sum([1 for x in clusteredSeqs if x[1] == -1])} sequences unclustered.")
  #logger.info("Writing down cluster information.")
  write_cluster(outdir, clusteredSeqs, inputSequences, goi)
  if goi:
    for header, cluster in goi2Cluster.items():
      logger.info(f"You'll find {header} (genome of interest) in cluster {cluster}.")

  logger.info("Extracting centroid sequences.")
  centroids = determine_centroids(clusteredSeqs, outdir, proc)

  logger.info("Generating clustered sequence file.")
  report_centroids(inputSequences, centroids, goi, outdir, clusteredSeqs)

  logger.info("Cleaning temporary files.")
  clean_up(outdir)

  sys.exit(0)

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def __create_logger():
    """
    doc string.
    """

    logger = logging.getLogger()
    #logger.setLevel(logging.WARNING)

    handle = logging.StreamHandler()
    #handle.setLevel(logging.WARNING)

    formatter = ColoredFormatter("%(log_color)sViralClust %(levelname)s -- %(asctime)s -- %(message)s", "%Y-%m-%d %H:%M:%S",
                                    log_colors={
                                            'DEBUG':    'bold_cyan',
                                            'INFO':     'bold_white',
                                            'WARNING':  'bold_yellow',
                                            'ERROR':    'bold_red',
                                            'CRITICAL': 'bold_red'}
                                )

    handle.setFormatter(formatter)
    logger.addHandler(handle)
    return logger

def create_outdir(outdir):
    try:
      os.makedirs(outdir)
      logger.info(f"Creating output directory: {outdir}")
    except FileExistsError:
      logger.warning(f"The output directory exists. Files will be overwritten.")

def __parse_arguments(d_args):
  """
  Parse all given arguments and check for error (e.g. file names).

  Arguments:
  d_args -- dict with input parameters and their values

  Returns:
  parsed and checked parameter list
  """

  if d_args['--version']:
    print("ViralClust version 0.1")
    exit(0)


  verbose = d_args['--verbose']
  if verbose:
    logger.setLevel(logging.INFO)

  inputSequences = d_args['<inputSequences>']
  if not os.path.isfile(inputSequences):
    logger.error("Couldn't find input sequences. Check your file")
    exit(1)

  goi = d_args['<genomeOfInterest>']
  if goi and not os.path.isfile(goi):
    logger.error("Couldn't find genome of interest. Check your file")
    exit(1)

  try:
    KMER = int(d_args['--kmer'])
  except ValueError:
    logger.error("Invalid parameter for k-mer size. Please input a number.")
    exit(2)

  try:
    proc = int(d_args['--process'])
  except ValueError:
    logger.error("Invalid number for CPU cores. Please input a number.")
    exit(2)

  output = d_args['--output']
  if output == 'pwd':
    output = os.getcwd()
  now = str(datetime.now()).split('.')[0].replace(' ','_').replace(':','-')
  create_outdir(output)

  METRICES = [
              'euclidean', 'manhattan', 'chebyshev',
              'minkwoski', 'canberra', 'braycurtis',
              'mahalanobis', 'wminkowski',
              'seuclidean', 'cosine'
  ]
  metric = d_args['--metric']
  if metric not in METRICES:
    log.warn(f"Invalid metric chosen. Will default back to cosine distance.")
    metric = 'cosine'
  
  try:
    neighbors = int(d_args['--neighbors'])
    if neighbors <= 0:
      raise ValueError
      
  except ValueError:
      log.error("Invalid parameter for --neighbors. Please input a positive integer.")
      exit(2)
  
  try:
    threshold = float(d_args['--dThreshold'])
    if not 0.0 <= threshold < 1:
      raise ValueError
  except ValueError:
    log.error("Invalid parameter for --dThreshold. Please input a number between [0,1).")
    exit(2)

  try:
    dimension = int(d_args['--dimension'])
    if dimension < 1:
      raise ValueError
  except ValueError:
      log.error("Invalid parameter for --dimension. Please input a positive integer.")
      exit(2)
  
  
  try:
    clusterSize = int(d_args['--clusterSize'])
    if clusterSize < 1:
      raise ValueError
  except ValueError:
      log.error("Invalid parameter for --clusterSize. Please input a positive integer.")
      exit(2)

  if d_args['--minSample'] == "CLUSTERSIZE":
    minSample = clusterSize
  else:
    try:
      minSample = int(d_args['--minSample'])
      if minSample < 1:
        raise ValueError
    except ValueError:
        log.error("Invalid parameter for --minSample. Please input a positive integer.")
        exit(2)

  umap_flag = d_args['--umap']

  return (inputSequences, goi, output,  KMER, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, umap_flag)

def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. No subcluster will be created.")
    del clusterObject

def __parse_fasta(filePath):
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

def __find_record_by_header(filePath, goi, headerToFind):
  """
  """
  for header, seq in __parse_fasta(filePath):
    if header == headerToFind:
      return(seq)
  
  if goi:
    for header, seq in __parse_fasta(goi):
      if header == headerToFind:
        return(seq)

def __profileA(sequence, allKmers, KMER):
  """
  """

  profile = np.zeros(shape=(len(allKmers)), dtype=np.float32)
  for k in iter([sequence[start : start + KMER] for start in range(len(sequence) - KMER)]):
      try:
        profile[allKmers[k]] += 1
      except KeyError:
        continue
  try:
    profile = profile / np.sum(profile)
    return(profile)
  except ZeroDivisionError:
    print(id2header[header] + " skipped, due to too many N's.")
    return(None)

def sequences_to_profiles(sequences, outdir, allKmers, KMER, previousID=-1):
  """
  """
  mode = 'wb' if previousID == -1 else 'ab'
  with open(f"{outdir}/profiles.pkl", mode) as outputStream:
    idHead = previousID
    for header, seq in __parse_fasta(sequences):
      idHead += 1
      id2header[idHead] = header
      header2id[header] = idHead
      if mode == 'ab':
        goiHeader.append(idHead)

      pickle.dump((idHead,__profileA(seq, allKmers, KMER)), outputStream)
      #outputStream.write(str(idHead)+"\t"+'\t'.join(map(str,__profileA(seq, allKmers, KMER)))+"\n")
  return(idHead if idHead >= 21 else 1)

def __load_profiles(outdir, header=[]):
  """
  """
  profiles = []
  with open(f"{outdir}/profiles.pkl", "rb") as inputStream:
    try:
      while True:
        data = pickle.load(inputStream)
        if header and not data[0] in header:
          # print(header, data[0])
          continue
        profiles.append(data)
    except EOFError:
      return(profiles)

def reduce_with_pca(vector):
  """
  """
  pca_model = PCA()
  pca_model.fit(vector)
  variances = pca_model.explained_variance_ratio_
  for i, var in enumerate(variances):
    if sum(variances[:i]) > 0.7 or i > 50:
      break
  pca_model = PCA(i)
  clusterable_embedding = pca_model.fit_transform(vector)  
  return(clusterable_embedding)

def reduce_with_umap(vector, neighbors, threshold, metric):
  """
  """
  clusterable_embedding = umap.UMAP(
          n_neighbors=neighbors,
          min_dist=threshold,
          n_components=50,
          random_state=42,
          metric=metric,
      ).fit_transform(vector)
  return(clusterable_embedding)

def dimension_reduction(outdir, neighbors, threshold, metric, umap_flag):
  """
  """
  profiles = __load_profiles(outdir)
  vector = [x[1] for x in profiles]

  
  if not umap_flag:
    embedding = reduce_with_pca(vector)
  else:
    embedding = reduce_with_umap(vector, neighbors, threshold, metric)

  del vector
  del profiles
  gc.collect()

  return(embedding)

def cluster_with_hdbscan(embedding, metric):
  """
  """

  clusterer = hdbscan.HDBSCAN(core_dist_n_jobs=1)
  if metric == "cosine":
    embedding = normalize(embedding,norm='l2')
    clusterer.fit(embedding)
  else:
    clusterer.fit(embedding, metric=metric)

  clusterlabel = clusterer.labels_
  probabilities = clusterer.probabilities_

  allCluster = list(zip([x for x in id2header.keys()], clusterlabel))
  for i in set(clusterlabel):
    for idx, label in allCluster:
      if label == i:
        if idx in goiHeader:
          goi2Cluster[id2header[idx]] = i
  return(allCluster)

def write_cluster(outdir, allCluster, inputSequences, goi): 
  """
  """
  
  clusterlabel = [x[1] for x in allCluster]
  #with open(f'{outdir}/cluster.txt', 'w') as outStream:
  for i in set(clusterlabel):
    if not i == -1: continue
    with open(f'{outdir}/cluster{i}.fasta', 'w') as fastaOut:
      #  outStream.write(f">Cluster {i}\n")
      for idx, label in allCluster:
        if label == i:
          if idx in goiHeader:
            goi2Cluster[id2header[idx]] = i
        #  outStream.write(f"{id2header[idx]}\n")
          sequence = __find_record_by_header(inputSequences, goi, id2header[idx])
          fastaOut.write(f">{id2header[idx]}\n{sequence.split('X'*10)[0]}\n")
      # outStream.write("\n")

def calc_pd(seqs):
  """
  """
  for element in seqs:
    try:
      stuff = (element[0], element[1])
    except TypeError:
      return None
  seq1, *profile1 = seqs[0]
  seq2, *profile2 = seqs[1]
  dFunction = scipyDistances[metric]
  distance = dFunction(profile1, profile2)
  return (seq1, seq2, distance)

def determine_centroids(allCluster, outdir, proc):
  """
  """

  multiPool = Pool(processes=proc)
  centroids = []
  seqCluster = { x : [] for x in set([x[1] for x in allCluster]) }


  for idx, cluster in allCluster:
    seqCluster[cluster].append(idx)

  

  for cluster, sequences in seqCluster.items():
    if cluster == -1:
      continue

    subProfiles = __load_profiles(outdir, sequences)
    indexForMatrix = {}
    seqHeader2Index = {}
    tmp = []
    for idx, profile in enumerate(subProfiles):
      indexForMatrix[idx] = profile[0]
      seqHeader2Index[profile[0]] = idx
      profile = (idx,) + profile[1:]
      tmp.append(profile)
    subProfiles = tmp[:]
    del tmp
      
    #indexForMatrix = { i : header for enumerate([x[0] for x in subProfiles]) }
    matrix = np.ones(shape=(len(indexForMatrix)+1,len(indexForMatrix)+1), dtype=np.float32)
    try:
      for result in multiPool.map(calc_pd, itertools.combinations(subProfiles, 2)):
        seq1, seq2, dist = result
        matrix[seq1][seq2] = dist
        matrix[seq2][seq1] = dist
    except IndexError:
      print(seq1, seq2, len(indexForMatrix),)
      print(len(matrix))
      exit(0)

    tmpMinimum = math.inf
    centroidOfCluster = -1

    if len(sequences) == 1:
      centroidOfCluster = cluster[0]
      centroids.append(centroidOfCluster)
      continue

    for sequence in sequences:
      averagedDistance = 0

      if sequence in goiHeader:
        continue

      for neighborSequence in sequences:
        if sequence == neighborSequence:
          continue
        averagedDistance += matrix[seqHeader2Index[sequence]][seqHeader2Index[neighborSequence]]

      averagedDistance /= len(sequences)-1

      if averagedDistance < tmpMinimum:
        tmpMinimum = averagedDistance
        centroidOfCluster = sequence

    centroids.append(centroidOfCluster)
  multiPool.close()
  multiPool.join()
  return(centroids)

def report_centroids(inputSequences, centroids, goi, outdir, allCluster):
  """
  """
  outputPath = f'{outdir}/{os.path.splitext(os.path.basename(inputSequences))[0]}_hdbscan.fasta'
  with open(outputPath, 'w') as outputStream:
    for centroidHeader in centroids:
      representativeSequence = __find_record_by_header(inputSequences, goi, id2header[centroidHeader])
      representativeSequence = representativeSequence.split('X'*10)[0]
      outputStream.write(f">{id2header[centroidHeader]}\n{representativeSequence}\n")
    # if goi:
    #   with open(goi, 'r') as goiStream:
    #     for line in goiStream:
    #       if line.startswith(">"):
    #         header = line.rstrip("\n").replace(':','_').replace(' ','_').strip(">_")
    #         if header2id[header] in centroids:
    #           continue
    #         outputStream.write(f">{header}\n{__find_record_by_header(inputSequences, goi, header)}\n")
  
  clusterlabel = [x[1] for x in allCluster]
  with open(f'{outputPath}.clstr', 'w') as outStream:
    for i in set(clusterlabel):
      outStream.write(f">Cluster {i}\n")
      seqInCluster = set([idx for idx,label in allCluster if label == i])
      for idx, seqIdx in enumerate(seqInCluster):
        sequence = __find_record_by_header(inputSequences, goi, id2header[seqIdx])
        outStream.write(f"{idx}\t{len(sequence)}nt, >{id2header[seqIdx]} ")
        if seqIdx in centroids:
          outStream.write('*\n')
        elif seqIdx in goiHeader:
          outStream.write("GOI\n")
        else:
          outStream.write("at +/13.37%\n")

def clean_up(outdir):
  """
  """
  os.remove(f"{outdir}/profiles.pkl")


if __name__ == "__main__":
  main()

