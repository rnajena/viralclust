#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

from Bio import SeqIO
import time
import sys
import pickle
import glob
import re


files = glob.glob(f"{sys.argv[1]}/gbvrl*.seq")
gb_vrl =  SeqIO.index_db(f"{sys.argv[1]}/gbvrl.idx", files, "genbank")


countryRegex = re.compile(r'country="([^"]+)"')
accessionDateRegex = re.compile(r'collection_date="([^"]+)"')

d_metaInformation = {}

for accession in gb_vrl:
  genbankEntry = gb_vrl.get_raw(accession).decode("utf-8")

  country = re.findall(countryRegex, genbankEntry)
  if country:
    country = country[0].replace('\n',' ').replace(' ','_').replace(',','')
  else:
    country = '--'

  accessionDate = re.findall(accessionDateRegex, genbankEntry)
  if accessionDate:
    accessionDate = accessionDate[0].replace('\n',' ').replace(' ','_').replace(',','')
  else:
    accessionDate = '--'

  taxonomy = gb_vrl[accession].annotations['taxonomy'] + [gb_vrl[accession].annotations['organism']]
  taxonomy = list(map(lambda x : x.replace('\n',' ').replace(' ','_'), taxonomy))

  family = 'unclassified family'

  for idx, element in enumerate(taxonomy):
    if element.endswith("viridae"):
      family = element
      familyIDX = idx
      break


  if taxonomy[idx+1:-1]:
    for element in taxonomy[idx+1:-1]:
      if element.endswith("virus"):
        genus = element
        break
  else:
    if family != "unclassified":
      genus = f'unclassified {family}'
    else:
      genus = "unclassified"

  organism = taxonomy[-1]
  taxonomy = family, genus, organism

  d_metaInformation[accession.split('.')[0]] = (country.strip(), accessionDate.strip(), taxonomy)

currentTime = time.asctime()
pickle.dump((currentTime, d_metaInformation), open(f'{sys.argv[1]}/ncbi_metainfo.pkl', 'wb'))