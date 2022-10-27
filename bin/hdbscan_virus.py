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

Other Dependencies:
  cd-hit

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
                                          'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                          'canberra', 'braycurtis',
                                          'mahalanobis', 'wminkowski', 'seuclidean',
                                          'cosine'.
                                          If an invalid metric is set, ViralClust will default back to 
                                          the cosine distance.
                                          [Default: cosine]

  --pca                                  Flag that determines whether (instead of UMAP) a PCA is used for dimension reduction. [Default: False]

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

import numpy as np
from multiprocessing import Pool

from datetime import datetime

from colorlog import ColoredFormatter
from docopt import docopt
from Bio import Phylo

from collections import Counter
import multiprocessing as mp
from multiprocessing import shared_memory
import itertools
import math
import subprocess

import scipy
import random

import umap.umap_ as umap
import hdbscan
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA


import linecache
import os
import tracemalloc

def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


def __parse_fasta(filePath, goi=False):
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

def read_sequences():
  """
  """
  global dim
  global matrix 
  global d_profiles
  global shm_name


  idHead = -1
  fastaContent = {}
  for header, sequence in __parse_fasta(inputSequences):
    #if not subCluster:
    idHead += 1
    id2header[idHead] = header
    header2id[header] = idHead
    fastaContent[idHead] = sequence
    #else:
      # fastaContent[header2id[header]] = sequence

  if genomeOfInterest:
    for header, sequence in __parse_fasta(genomeOfInterest):
        if header in header2id:
          idHead = header2id[header]
        else:
          idHead += 1
        goiHeader.append(idHead)
        id2header[idHead] = header
        header2id[header] = idHead
        fastaContent[idHead] = sequence

  # if not subCluster:
  dim = len(fastaContent)
  matrix = np.ones(shape=(dim, dim), dtype=float)

  #shm = shared_memory.SharedMemory(create=True, size=matrix.nbytes)
  #matrix = np.ndarray((dim, dim), dtype=float, buffer=shm.buf)

  if len(fastaContent) < 21:
    return 1
  return fastaContent


def profileA(entry):
  """
  """
  # global shm_name
  header, sequence = entry
  # existing_shm = shared_memory.SharedMemory(name=shm_name)
  # allProfiles = np.ndarray((dim, len(allKmers)), buffer=existing_shm.buf)
  # profile = allProfiles[header]
  profile = [0]*len(allKmers)
  #profile = np.empty(shape=(len(allKmers)), dtype=np.float64)
  for k in iter([sequence[start : start + KMER] for start in range(len(sequence) - KMER)]):
      try:
        profile[allKmers[k]] += 1
      except KeyError:
        continue
  kmerSum = np.sum(profile)
  try:
    profile = profile / kmerSum
    return (header, profile)
  except ZeroDivisionError:
    print(id2header[header] + " skipped, due to too many N's.")
    return(None, None)
  # exisiting_shm.close()
    #print(header, id2header[header], kmerSum)
    #print(sequence)
    #exit(1)

def determine_profile(proc):

  global d_sequences
  global d_profiles
  p = proc
  tracemalloc.start()
  for entry in d_sequences.items():
    header, profile = profileA(entry)
  #allProfiles = p.map(profileA, d_sequences.items())
  #p.close()
  #p.join()
  
  
  #for header, profile in allProfiles:
    if header:
      d_profiles[header] = profile    
  snapshot = tracemalloc.take_snapshot()
  display_top(snapshot)
  
def calc_pd(seqs):
  """
  """
  for element in seqs:
    try:
      stuff = (element[0], element[1])
    except TypeError:
      return None
  seq1, profile1 = seqs[0]
  seq2, profile2 = seqs[1]
  dFunction = scipyDistances[metric]
  distance = dFunction(profile1, profile2)
  return (seq1, seq2, distance)
  #

def apply_umap():
  """
  """
  global clusterlabel
  global probabilities
  global allCluster
  global d_profiles
  
  profiles = [(idx,profile) for idx, profile in d_profiles.items() if idx in d_sequences]
  vector = [x[1] for x in profiles]

  try:
    
    if pca_flag:
      pca_model = PCA()
      pca_model.fit(vector)
      variances = pca_model.explained_variance_ratio_
      for i, var in enumerate(variances):
        if sum(variances[:i]) > 0.7 or i > 50:
          break
      #for i in range(len(pca_model.explained_variance_)):
      #  if sum(pca_model.explained_variance_ratio_)
      pca_model = PCA(i)
      clusterable_embedding = pca_model.fit_transform(vector)        
    else:
      clusterable_embedding = umap.UMAP(
            n_neighbors=neighbors,
            min_dist=threshold,
            n_components=50,
            random_state=42,
            metric=metric,
        ).fit_transform(vector)

    clusterer = hdbscan.HDBSCAN()
    if metric == "cosine":
      clusterable_embedding = normalize(clusterable_embedding,norm='l2')
      clusterer.fit(clusterable_embedding)
    else:
      clusterer.fit(clusterable_embedding,metric=metric)

    clusterlabel = clusterer.labels_
    probabilities = clusterer.probabilities_
    if len(set(clusterlabel)) == 1:
      raise TypeError

  except TypeError:
    import shutil
    shutil.copyfile(inputSequences, f'{outdir}/{os.path.splitext(os.path.basename(inputSequences))[0]}_hdbscan.fasta')
    return 1

  allCluster = list(zip([x[0] for x in profiles], clusterlabel))
  #if not subCluster:
  with open(f'{outdir}/cluster.txt', 'w') as outStream:
    for i in set(clusterlabel):
      with open(f'{outdir}/cluster{i}.fasta', 'w') as fastaOut:
        outStream.write(f">Cluster {i}\n")
        for idx, label in allCluster:
          if label == i:
            if idx in goiHeader:
              goi2Cluster[id2header[idx]] = i
            outStream.write(f"{id2header[idx]}\n")
            fastaOut.write(f">{id2header[idx]}\n{d_sequences[idx].split('X'*10)[0]}\n")
      outStream.write("\n")

def get_centroids(proc):
  """
  """

  global clusterlabel
  global allCluster
  global d_profiles
  global matrix

  centroids = []

  seqCluster = { x : [] for x in set(clusterlabel)}

  for idx, cluster in allCluster:
    seqCluster[cluster].append(idx)

  p = proc
  
  for cluster, sequences in seqCluster.items():
    if cluster == -1:
      #for sequence in sequences:
        #if sequence in goiHeader:
        #  centroids.append(sequence)
      continue #print(sequence)

    subProfiles = {seq : profile for seq, profile in d_profiles.items() if seq in sequences}

    #if not subCluster:
    for result in p.map(calc_pd, itertools.combinations(subProfiles.items(), 2)):
      seq1, seq2, dist = result
      matrix[seq1][seq2] = dist
      matrix[seq2][seq1] = dist

    tmpMinimum = math.inf
    centroidOfCluster = -1

    if len(sequences) == 1:
      centroidOfCluster = cluster[0]
      centroids.append(centroidOfCluster)
      continue

    for sequence in sequences:
      averagedDistance = 0

      if sequence in goiHeader:
        #print(goiHeader)
        #centroids.append(sequence)
        continue

      for neighborSequence in sequences:
        if sequence == neighborSequence:
          continue
        averagedDistance += matrix[sequence][neighborSequence]

      averagedDistance /= len(sequences)-1

      if averagedDistance < tmpMinimum:
        tmpMinimum = averagedDistance
        centroidOfCluster = sequence

    centroids.append(centroidOfCluster)
  p.close()
  p.join()
  return(centroids)

def output_centroids(centroids):
  """
  """

  reprSeqs = {centroid : d_sequences[centroid] for centroid in centroids}
  reprSeqs.update({goiID : d_sequences[goiID] for goiID in goiHeader})
  outputPath = f'{outdir}/{os.path.splitext(os.path.basename(inputSequences))[0]}_hdbscan.fasta'

  with open(outputPath, 'w') as outStream:
    for centroidID, sequence in reprSeqs.items():
      outStream.write(f">{id2header[centroidID]}\n{sequence}\n")

  with open(f'{outputPath}.clstr', 'w') as outStream:
    for i in set(clusterlabel):
      outStream.write(f">Cluster {i}\n")
      seqInCluster = set([idx for idx,label in allCluster if label == i])
      for idx, seqIdx in enumerate(seqInCluster):
        outStream.write(f"{idx}\t{len(d_sequences[seqIdx])}nt, >{id2header[seqIdx]} ")
        if seqIdx in centroids:
          outStream.write('*\n')
        elif seqIdx in goiHeader:
          outStream.write("GOI\n")
        else:
          outStream.write("at +/13.37%\n")

###################################################################################################

inputSequences = None
goi = None
outdir = None
k = None
proc = None
OUT = ''

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn


def create_logger():
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
      #os.makedirs(f"{outdir}/tmpSequences")
      logger.info(f"Creating output directory: {outdir}")
    except FileExistsError:
      logger.warning(f"The output directory exists. Files will be overwritten.")

def parse_arguments(d_args):
  """
  Parse all given arguments and check for error (e.g. file names).

  Arguments:
  d_args -- dict with input parameters and their values

  Returns:
  parsed and checked parameter list
  """

  if d_args['--version']:
    print("viralClust version 0.1")
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
  #output = f"{output}/viralClust-{now}"
  create_outdir(output)

  METRICES = [
              'euclidean', 'manhatten', 'chebyshev',
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

  pca_flag = d_args['--pca']

  return (inputSequences, goi, output,  KMER, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, pca_flag)

def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. No subcluster will be created.")
    del clusterObject

def perform_clustering():
  
  global clusterlabel
  global d_sequences
  global shm_name

  d_sequences = read_sequences()
  if d_sequences == 1:
    import shutil
    shutil.copyfile(inputSequences, f'{outdir}/{os.path.splitext(os.path.basename(inputSequences))[0]}_hdbscan.fasta')
    return 1

  d_profiles = {}
  # d_profiles = np.ones(shape=(dim,len(allKmers)), dtype=np.float64)
  # shm = shared_memory.SharedMemory(create=True, size=d_profiles.nbytes)
  # shm_name = shm.name
  # logger.warn(shm_name)
  # d_profiles = np.ndarray(d_profiles.shape, dtype=np.float64, buffer=shm.buf)

  multiPool = Pool(processes=proc)
  #virusClusterer = Clusterer(logger, inputSequences, outdir,  k, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, pca_flag, goi=goi)


  logger.info("Determining k-mer profiles for all sequences.")
  determine_profile(multiPool)
  if goi:
    logger.info(f"Found {len(virusgoiHeader)} genome(s) of interest.")
  logger.info("Clustering with UMAP and HDBSCAN.")
  code = apply_umap()
  if code == 1:
    __abort_cluster(inputSequences)
    return 0
  clusterInfo = clusterlabel
  #logger.info(f"Summarized {dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")

  goiCluster = goi2Cluster
  if goiCluster:
    for header, cluster in goiCluster.items():
      logger.info(f"You find the genome {header} in cluster {cluster}.")

  logger.info("Extracting centroid sequences and writing results to file.\n")
  multiPool = Pool(processes=proc)
  centroids = get_centroids(multiPool)
  output_centroids(centroids)

  logger.info(f"All done.")
  # sequences = d_sequences
  # distanceMatrix = matrix
  # profiles = d_profiles
  #virusshm.close()
  #del virusClusterer
  # shm.close()
  # shm.unlink()
  return 0


if __name__ == "__main__":
  
  logger = create_logger()
  (inputSequences, goi, outdir,  KMER, proc, metric, neighbors, threshold, dimension, clusterSize, minSample, pca_flag) = parse_arguments(docopt(__doc__))

  reducedSequences = f"{outdir}/{os.path.basename(inputSequences)}"
  nucleotides = set(["A","C","G","T"])
  allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(nucleotides, repeat=KMER))}

  id2header = {}
  header2id = {}
  dim = 0
  matrix = np.empty(shape=(dim,dim))

  genomeOfInterest = ''
  goiHeader = []
  goi2Cluster = {}
  shm_name = ""
  d_profiles = {}

  scipyDistances = {
      'euclidean' : scipy.spatial.distance.euclidean ,
      'manhatten' : scipy.spatial.distance.cityblock ,
      'chebyshev' : scipy.spatial.distance.chebyshev  ,
      'minkwoski': scipy.spatial.distance.minkowski ,
      'canberra' : scipy.spatial.distance.canberra ,
      'braycurtis' : scipy.spatial.distance.braycurtis ,
      'mahalanobis' : scipy.spatial.distance.mahalanobis ,
      'wminkowski' : scipy.spatial.distance.wminkowski ,
      'seuclidean' : scipy.spatial.distance.seuclidean ,
      'cosine' : scipy.spatial.distance.cosine  
  }

  #d_sequences = {}
  #centroids = []
  #allCluster = []
  #clusterlabel = []
  #probabilities = []
  
  if goi:
    genomeOfInterest = goi

  logger.info("Starting to cluster you data. Stay tuned.")
  perform_clustering()
  logger.info("Linking the latest results.")
  if os.path.islink(f"{os.path.dirname(outdir)}/latest"):
    os.remove(f"{os.path.dirname(outdir)}/latest")
  os.system(f"ln -s {outdir} {os.path.dirname(outdir)}/latest")
