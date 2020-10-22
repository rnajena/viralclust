#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

import re
from Bio import Entrez
from collections import defaultdict
import time
import sys
import pickle

import utils

sequenceFile = sys.argv[1]

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
  accID2Name = {accessionID : scientificNames[taxID].replace(' ','_') if taxID in scientificNames else "unclassified" for accessionID, taxID in taxIDs.items()}
  return accID2Name

def retrieve_taxonomy(prefix):
  Entrez.email = "some_example@mail.com"
  genbankACCRegex = re.compile(r'([A-Z]{2}[0-9]{8}|[A-Z]{2}[0-9]{6}|[A-Z][0-9]{5}|NC_[0-9]{6})')
  accessionIDs = [header for header, sequence in utils.parse_fasta(sequenceFile)]

  validIDs = list(filter(genbankACCRegex.match, accessionIDs))
  invalidIDs = [x for x in accessionIDs if x not in validIDs]

  while 1:
    try:
      handle = Entrez.elink(dbfrom='nuccore', db='taxonomy', id=validIDs, idtype='acc', rettype='xml')
      record = Entrez.read(handle)
      handle.close()
    except RuntimeError:
      time.sleep(5)
      continue
    break


  taxIDs = { str(x['IdList'][0]).split('.')[0] : x['LinkSetDb'][0]['Link'][0]['Id'] if x['LinkSetDb'] else "XXXXX" for x in record }
  values = list(set(taxIDs.values()))

  target_handle = Entrez.efetch(db='taxonomy', id=','.join(values), retmode='xml')
  target_record = Entrez.read(target_handle)
  target_handle.close()

  accID2species = get_clade_by_rank('species', taxIDs, target_record)
  accID2genus = get_clade_by_rank('genus', taxIDs, target_record)

  accID2GBdescription = {}
  nuccore_handle = Entrez.efetch(db='nuccore', id=','.join(list(accID2species.keys())), rettype='gb', retmode='xml')
  nuccore_records = Entrez.read(nuccore_handle)

  for record in nuccore_records:
    country = collectionDate = '--'
    accID = record['GBSeq_primary-accession'].split('.')[0]
    featureTable = record['GBSeq_feature-table'][0]['GBFeature_quals']
    for feature in featureTable:
      if feature['GBQualifier_name'] == 'country':
        country = feature['GBQualifier_value']
      if feature['GBQualifier_name'] == 'collection_date':
        collectionDate = feature['GBQualifier_value']
    accID2GBdescription[accID] = f"{accID2species[accID]},{accID2genus[accID]},{record['GBSeq_definition'].replace(' ','_')},{country.replace(' ','_')},{collectionDate.replace(' ','_')}"

  while len(accID2GBdescription) != len(accID2species):
    remainder = [acc for acc in accID2species if acc not in accID2GBdescription]
    time.sleep(1)
    nuccore_handle = Entrez.efetch(db='nuccore', id=','.join(remainder), rettype='gb', retmode='xml')
    nuccore_records = Entrez.read(nuccore_handle)
    for record in nuccore_records:
      country = collectionDate = '--'
      accID = record['GBSeq_primary-accession'].split('.')[0]
      featureTable = record['GBSeq_feature-table'][0]['GBFeature_quals']
      for feature in featureTable:
        if feature['GBQualifier_name'] == 'country':
          country = feature['GBQualifier_value']
        if feature['GBQualifier_name'] == 'collection_date':
          collectionDate = feature['GBQualifier_value']
      accID2GBdescription[accID] = f"{accID2species[accID]},{accID2genus[accID]},{record['GBSeq_definition'].replace(' ','_')},{country.replace(' ','_')},{collectionDate.replace(' ','_')}"
  return accID2GBdescription

accID2Desc = retrieve_taxonomy(sequenceFile)

with open('ncbiMETA.pkl', 'wb') as outputStream:
  pickle.dump(accID2Desc, outputStream)
