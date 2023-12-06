#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

# ToDo: ViralClust Summary csv
# ToDo: Optional -- NCBI Taxonomy Info
# ToDo: Clustered/UnClustered Sequences linken
# ToDo: Herausfinden, wie JavaScript funktioniert: RadioButtons fuer GOIs


# * Default imports
import sys
import os
import glob
import json
import colorsys
import random
from collections import defaultdict

# * 3rd Party imports, need to be installed
import numpy as np
import seaborn as sns

# * Own file imports, should be in the same directory
import html_parser

HTMLHEADER = ''
HTMLBODY = ''


def main():
  """

  """


  # fastaFile = sys.argv[1]
  clusterInfoDir = sys.argv[1]
  HTMLREPORT = f"{sys.argv[2]}viralclust.html"
  global HTMLBODY

  prepare_html_document()

  tools = ["cdhitest", "sumaclust", "vclust", "mmseqs", "hdbscan"]
  toolClusterInformation, toolCentroids, toolUnclustered, acc2toolCluster = fill_informations(tools, clusterInfoDir)

  maxNumberCluster = max([len(x) for x in toolCentroids.values()])
  snsPalettes = ["pastel", "muted", "deep", "bright", "dark"]
  ratio = int(maxNumberCluster/len(snsPalettes))+1
  sampledColors = [x for p in snsPalettes for x in sns.color_palette(p, ratio).as_hex() ]
  
  id2color = { clusterID : color for clusterID, color in zip(range(0,maxNumberCluster), sampledColors) }
  id2color.update({-1 : 'ffffff'})

  for accessionID, tool2clusterID in acc2toolCluster.items():
    HTMLBODY += f'<tr>\n<td align="center"><a href="https://www.ncbi.nlm.nih.gov/nuccore/{accessionID}">{accessionID}</a></td>\n'
    for tool in tools:
      try:
        if tool2clusterID[tool] == -1:
          raise KeyError
        HTMLBODY += f'<td align="center", border="2", bgcolor={id2color[tool2clusterID[tool]]}>Cluster {tool2clusterID[tool]}</td>\n'
      except KeyError:
        HTMLBODY += f'<td align="center">Unclustered</td>\n'
    HTMLBODY += f'</tr>\n'
      

  close_html_body()
  print_html_to_file(HTMLREPORT)
  sys.exit(0)  

def prepare_html_document():
  """

  """
  global HTMLHEADER
  global HTMLBODY
  HTMLHEADER = """
  <!doctype html>
  <html lang="en">

  <head>
      <meta http-equiv="Content-type" content="text/html; charset=utf-8">
      <meta name="viewport" content="width=device-width,initial-scale=1,user-scalable=no">
      <title>ViralClust</title>
      <link rel="stylesheet" type="text/css" href="jquery.dataTables.min.css">

      <style type="text/css" media="screen">
          @import url('//cdn.datatables.net/1.10.2/css/jquery.dataTables.css');
          td.details-control {
              background: url('http://www.datatables.net/examples/resources/details_open.png') no-repeat center center;
              cursor: pointer;
          }

          tr.shown td.details-control {
              background: url('http://www.datatables.net/examples/resources/details_close.png') no-repeat center center;
          }
      </style>

      <script type="text/javascript" language="javascript" src="jquery-3.3.1.min.js"></script>
      <script type="text/javascript" language="javascript" src="jquery.dataTables.min.js"></script>
      <script type="text/javascript" class="init">

          function format(value) {
              return '<tr>' + value + '</tr>';
          }
          $(document).ready(function () {
              var table = $('#viralclust').DataTable({
                paging: false,
                orderClasses: false
              });
          });


      </script>
  </head>
  """
              
  HTMLBODY = """
  <body>

      <table id="viralclust" class="display nowrap" cellspacing="0" width="50%", border="1px">
          <thead>
              <tr>
                  <th>Accession ID</th>    
                  <th>cd-hit-est</th>
                  <th>sumaclust</th>
                  <th>vclust</th>
                  <th>MMseqs2</th>
                  <th>HDBSCAN</th>
              </tr>
          </thead>
    <tbody>
  """

def close_html_body():
  """

  """
  global HTMLBODY

  HTMLBODY += """
        </tbody>
    </table>

    </body>
    </html>
  """

def print_html_to_file(HTMLREPORT):
  """
  """
  with open(HTMLREPORT, 'w') as outputStream:
    outputStream.write(HTMLHEADER + "\n")
    outputStream.write(HTMLBODY + "\n")

def fill_informations(tools, clusterInfoDir):
  """

  """

  toolClusterInformation = defaultdict(list)
  toolCentroids = defaultdict(list)
  toolUnclustered = defaultdict(list)
  acc2toolCluster = defaultdict(dict)

  for tool in tools:
    clusterInfoFile = glob.glob(f"{clusterInfoDir}/*{tool}*.clstr")[0]
    acc2cluster, centroids, failbob, cluster, gois = html_parser.parse_clusterFile(clusterInfoFile)
    toolClusterInformation[tool] = acc2cluster
    toolCentroids[tool] = centroids
    toolUnclustered[tool] = failbob
    for accID, clusterID in acc2cluster.items():
      acc2toolCluster[accID][tool] = clusterID
  
  return(toolClusterInformation, toolCentroids, toolUnclustered, acc2toolCluster)

if __name__ == "__main__":
  main()