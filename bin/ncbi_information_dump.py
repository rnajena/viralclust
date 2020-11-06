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
    country = country[0]
  else:
    country = '--'

  accessionDate = re.findall(accessionDateRegex, genbankEntry)
  if accessionDate:
    accessionDate = accessionDate[0]
  else:
    accessionDate = '--'

  taxonomy = gb_vrl[accession].annotations['taxonomy'][-2:]


  d_metaInformation[accession.split('.')[0]] = (country, accessionDate, taxonomy)

currentTime = time.asctime()

pickle.dump((currentTime, d_metaInformation), open('test.pkl', 'wb'))