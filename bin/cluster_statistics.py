#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
Usage:
  cluster_statistics.py [options] <treeFile> <seqFile> <clusterFile>

Options:
  -h, --help        Prints this help and exits.
  --ncbi            Connects to the NCBI nuccore and taxonomy database to get taxonomic information based of the
                    header line of each fasta record. Note: Recommended, if a GenBank ID is part of the header.
  --toolName NAME   Used as prefix for NCBI taxonomy summary.

"""

import sys
import dendropy
import re
import numpy as np
import time

import utils

from docopt import docopt
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict

def get_clade_by_rank(rank, taxIDs, target_record):
  scientificNames = {}

  for x in target_record:
    if x['Rank'] == rank:
      scientificNames[x['TaxId']] = x['ScientificName']
    else:
      for level in x['LineageEx']:
        if level['Rank'] == rank:
          scientificNames[x['TaxId']] = level['ScientificName']
          break
  accID2Name = {accessionID : scientificNames[taxID] if not taxID == 'XXXXX' else "unclassified" for accessionID, taxID in taxIDs.items()}
  return accID2Name

def retrieve_taxonomy(prefix):
  speciesPerCluster = 0
  genusPerCluster = 0
  clusterPerSpecies = defaultdict(set)
  clusterPerGenus = defaultdict(set)

  Entrez.email = "some_example@mail.com"
  genbankACCRegex = re.compile(r'([A-Z]{2}[0-9]{8}|[A-Z]{2}[0-9]{6}|[A-Z][0-9]{5}|NC_[0-9]{6})')

  with open(f'{prefix}_taxonomy_info.txt', 'w') as outputStream:
    for clusterID, accessionIDs in cluster.items():
      validIDs = list(filter(genbankACCRegex.match, accessionIDs))
      invalidIDs = [x for x in accessionIDs if x not in validIDs]

      handle = Entrez.elink(dbfrom='nuccore', db='taxonomy', id=validIDs, idtype='acc', rettype='xml')
      record = Entrez.read(handle)
      handle.close()

      taxIDs = { str(x['IdList'][0]).split('.')[0] : x['LinkSetDb'][0]['Link'][0]['Id'] if x['LinkSetDb'] else "XXXXX" for x in record }
      target_handle = Entrez.efetch(db='taxonomy', id=list(taxIDs.values()), retmode='xml')
      target_record = Entrez.read(target_handle)
      target_handle.close()

      accID2species = get_clade_by_rank('species', taxIDs, target_record)
      accID2genus = get_clade_by_rank('genus', taxIDs, target_record)

      if len(accessionIDs) > 1:
        speciesPerCluster += len(set(accID2species.values()))
        genusPerCluster += len(set(accID2genus.values()))
        for accID, species in accID2species.items():
          clusterPerSpecies[species].add(clusterID)
        for accID, genus in accID2genus.items():
          clusterPerGenus[genus].add(clusterID)

      accID2GBdescription = {}
      nuccore_handle = Entrez.efetch(db='nuccore', id=','.join(list(accID2species.keys())), rettype='gb', retmode='xml')
      nuccore_records = Entrez.read(nuccore_handle)
      for record in nuccore_records:
        country = collectionDate = '--'
        featureTable = record['GBSeq_feature-table'][0]['GBFeature_quals']
        for feature in featureTable:
          if feature['GBQualifier_name'] == 'country':
            country = feature['GBQualifier_value']
          if feature['GBQualifier_name'] == 'collection_date':
            collectionDate = feature['GBQualifier_value']
        # accID2GBdescription[gb_record.id.split('.')[0]] = gb_record.description
        accID2GBdescription[record['GBSeq_primary-accession'].split('.')[0]] = record['GBSeq_definition'] + ', ' + country + ', ' + collectionDate

      while len(accID2GBdescription) != len(accID2species):
        remainder = [acc for acc in accID2species if acc not in accID2GBdescription]
        time.sleep(1)
        nuccore_handle = Entrez.efetch(db='nuccore', id=','.join(remainder), rettype='gb', retmode='xml')
        nuccore_records = Entrez.read(nuccore_handle)
        for record in nuccore_records:
          country = collectionDate = '--'
          featureTable = record['GBSeq_feature-table'][0]['GBFeature_quals']
          for feature in featureTable:
            if feature['GBQualifier_name'] == 'country':
              country = feature['GBQualifier_value']

            if feature['GBQualifier_name'] == 'collection_date':
              collectionDate = feature['GBQualifier_value']

          accID2GBdescription[record['GBSeq_primary-accession'].split('.')[0]] = f"{record['GBSeq_definition'].replace(' ','_')},{country.replace(' ','_')},{collectionDate.replace(' ','_')}"
        #for gb_record in SeqIO.parse(nuccore_handle, 'genbank'):
        #  accID2GBdescription[gb_record.id.split('.')[0]] = gb_record.description
      # print(list(accID2species.keys()))
      # print(list(accID2GBdescription.keys()))
      # exit(0)
      #for accID in accID2species:
      #  if accID not in accID2GBdescription:
      #    accID2GBdescription[accID] = "No GenBank entry found"



      outputStream.write(f"Cluster: {int(clusterID)+1}\n")
      for acc, genusName in accID2genus.items():
        outputStream.write(f"{acc},{genusName.replace(' ','_')},{accID2species[acc].replace(' ','_')},{accID2GBdescription[acc]}\n")
      for acc in invalidIDs:
        outputStream.write(f"{acc}, --, --, --\n")
      outputStream.write(f"####################\n")



  avgClusterPerSpecies = sum([len(x) for x in clusterPerSpecies.values()])/len(clusterPerSpecies)
  avgClusterPerGenus = sum([len(x) for x in clusterPerGenus.values()])/len(clusterPerGenus)
  return(avgClusterPerSpecies, avgClusterPerGenus)

#########################################################################

args = docopt(__doc__)


treeFile = args['<treeFile>']
seqFile = args['<seqFile>']
clusterFile = args['<clusterFile>']
NCBI = args['--ncbi']
PREFIX = args['--toolName']

allSequences = {header : seq for header,seq in utils.parse_fasta(seqFile)}
cluster, centroids, failbob = utils.parse_clusterFile(clusterFile)

# HDBScan specific. -1 describe sequences that are not clustered at all.
#try:
#  cluster.pop('-1')
#except KeyError:
#  pass



if not failbob:
  failbob = [cluster for idx,cluster in cluster.items() if len(cluster) == 1]
realCluster = {idx : cluster for idx,cluster in cluster.items() if len(cluster) != 1 or idx != '-1'}
  # print(failbob)

tree = dendropy.Tree.get(path=treeFile, schema='newick')
dm = tree.phylogenetic_distance_matrix()

#print(tree.leaf_nodes())

overallSum = 0
for idx, cl in realCluster.items():
  clusterSum = 0
  for header in cl:
    if header in centroids:
      centroid = header
      centroid_taxon = tree.find_node_with_taxon_label(centroid.replace('_',' ')).taxon
      break
  for header in cl:
    if header == centroid:
      continue
    t1 = tree.find_node_with_taxon_label(header.replace('_',' ')).taxon
    clusterSum += (dm.distance(t1,centroid_taxon))
  overallSum += clusterSum / len(cl)
avgOverallSum = overallSum / len(realCluster)

allCluster = np.array([len(cl) for _,cl in realCluster.items()])

if NCBI:
  (avgClusterPerSpecies, avgClusterPerGenus) = retrieve_taxonomy(PREFIX)
else:
  avgClusterPerSpecies = avgClusterPerGenus = '--'




print(f"{len(allSequences)},{len(realCluster)},{np.min(allCluster)},{np.max(allCluster)},{np.mean(allCluster)},{np.median(allCluster)},{avgOverallSum},{len(failbob)},{avgClusterPerSpecies},{avgClusterPerGenus}")

