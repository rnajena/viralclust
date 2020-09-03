#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


import sys
import re

sequenceFile = sys.argv[1]

regex_orf = re.compile(r'M[^*]{25,}?\*')
genbankACCRegex = re.compile(r'[^A-Z]([A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{2}[0-9]{8}|NC_[0-9]{6})[^0-9]')

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


def rev_comp(sequence):
  """
  """
  d_comp = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
  return ("".join([d_comp[nt] if nt in d_comp else nt for nt in sequence[::-1]]))



def parse_fasta(filePath):
  """
  """

  with open(filePath, 'r') as inputStream:
    header = ''
    seq = ''

    for line in inputStream:
      if line.startswith(">"):
        if header:
          seq = seq + "X"*10 + rev_comp(seq)
          yield (header, seq)

        header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
        accessionID = set(re.findall(genbankACCRegex, header))
        if len(accessionID) == 1:
          header = accessionID.pop()
        seq = ''
      else:
        seq += line.rstrip("\n").upper().replace('U','T')

    seq = seq + "X"*10 + rev_comp(seq)
    yield (header, seq)

for header, sequence in parse_fasta(sequenceFile):
  positiveStrand = ""
  longestCDS = 0
  strands = sequence.split("X"*10)
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
  
  print(f">{header}\n{positiveStrand}")
  