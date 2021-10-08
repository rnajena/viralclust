#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


import sys
import re

import utils

sequenceFile = sys.argv[1]

regex_orf = re.compile(r'M[^*]{25,}?\*')

codon2aminoacid = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
  }



for header, sequence in utils.parse_fasta(sequenceFile):
  positiveStrand = sequence
  longestCDS = 0
  strands = [sequence, utils.reverseComplement(sequence)]
  for strand in strands:
    for frame in range(3):
      proteinSequence = ""
      for fragment in range(frame, len(strand), 3):
        codon = strand[fragment:fragment+3]
        if len(codon) != 3:
          continue
        try:
          proteinSequence += codon2aminoacid[codon]
        except KeyError:
          proteinSequence += 'X'
      matches = regex_orf.findall(proteinSequence)
      allORFs = "".join([x for x in matches if x])
      if len(allORFs)/float(len(strand)) > longestCDS:
        longestCDS = len(allORFs)/float(len(strand))
        positiveStrand = strand
  if positiveStrand:
    print(f">{header}\n{positiveStrand}")
