#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
viralClust v 0.1

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
  viralClust.py [options] <inputSequences> [<genomeOfInterest>]

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
  --version                               Prints the version of viralClust and exits.
  -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]

  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

  --subcluster                            Additionally to the initial cluster step, each cluster gets analyzed for 
                                          local structures and relations which results in subcluster for each cluster. [Default: False]

Version:
  viralClust v0.1 (alpha)
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

from ClusterViruses import Clusterer


def create_logger():
    """
    doc string.
    """

    logger = logging.getLogger()
    #logger.setLevel(logging.WARNING)

    handle = logging.StreamHandler()
    #handle.setLevel(logging.WARNING)

    formatter = ColoredFormatter("%(log_color)sVeGETA %(levelname)s -- %(asctime)s -- %(message)s", "%Y-%m-%d %H:%M:%S", 
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
    k = int(d_args['--kmer'])
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
  output = f"{output}/viralClust-{now}"
  create_outdir(output)

  subcluster = d_args['--subcluster']

  return (inputSequences, goi, output,  k, proc, subcluster)

def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. No subcluster will be created.")
    del clusterObject

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(logger, inputSequences, k, proc, outdir, goi=goi)

  logger.info("Removing 100% identical sequences.")
  code = virusClusterer.remove_redundancy()
  logger.info("Sequences are all parsed.")

  if code == 1:
    __abort_cluster(virusClusterer, inputSequences)
    return 0

  if goi:
    logger.info(f"Found {len(virusClusterer.goiHeader)} genome(s) of interest.")
  logger.info("Determining k-mer profiles for all sequences.")
  virusClusterer.determine_profile(multiPool)
  logger.info("Clustering with UMAP and HDBSCAN.")
  code = virusClusterer.apply_umap()
  if code == 1:
    __abort_cluster(virusClusterer, inputSequences)
    #logger.warning(f"All sequences fall into one cluster. Aligning this one without dividing the sequences anymore.")
    return 0
  clusterInfo = virusClusterer.clusterlabel
  logger.info(f"Summarized {virusClusterer.dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")

  goiCluster = virusClusterer.goi2Cluster
  if goiCluster:
    for header, cluster in goiCluster.items():
      logger.info(f"You find the genome {header} in cluster {cluster}.")

  logger.info("Extracting centroid sequences and writing results to file.\n")
  virusClusterer.get_centroids(multiPool)
  virusClusterer.split_centroids()
  
  logger.info(f"Extracting representative sequences for each cluster.")
  sequences = virusClusterer.d_sequences
  distanceMatrix = virusClusterer.matrix
  profiles = virusClusterer.d_profiles
  del virusClusterer

  if not subcluster:
    return 0

  for file in glob.glob(f"{outdir}/cluster*.fa"):
    if file == f"{outdir.rstrip('/')}/cluster-1.fa":
      continue
    virusSubClusterer = Clusterer(logger, file, k, proc, outdir, subCluster=True)
    code = virusSubClusterer.remove_redundancy()
    
    if code == 1:
      __abort_cluster(virusSubClusterer, file)
      #logger.warn(f"Too few sequences for clustering in {os.path.basename(file)}. Alignment will be calculated with all sequences of this cluster.")
      #del virusSubClusterer
      continue

    code = virusSubClusterer.apply_umap()
    
    if code == 1:
      __abort_cluster(virusSubClusterer, file)
      #logger.warn(f"Too few sequences for clustering in {os.path.basename(file)}. Alignment will be calculated with all sequences of this cluster.")
      #del virusSubClusterer
      continue

    virusSubClusterer.get_centroids(multiPool)
    virusSubClusterer.split_centroids()
    del virusSubClusterer

if __name__ == "__main__":
  logger = create_logger()
  (inputSequences, goi, outdir, k, proc, subcluster) = parse_arguments(docopt(__doc__))

  logger.info("Starting to cluster you data. Stay tuned.")
  perform_clustering()
  if os.path.islink(f"{os.path.dirname(outdir)}/latest"):
    os.remove(f"{os.path.dirname(outdir)}/latest")
  os.system(f"ln -s {outdir} {os.path.dirname(outdir)}/latest")

