
## ViralClust - Find representative viruses for your dataset
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.8](https://img.shields.io/badge/Language-Python_3.8-darkred.svg)![NextFlow](https://img.shields.io/badge/Nextflow-20.07.01-blue.svg)![conda](https://img.shields.io/badge/Uses-conda-green.svg)

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)

***

### DISCLAIMER
This pipeline is work-in-progress.
There are some bugs known to me, some aren't. Before getting desperate, please check out the Issues that are already opened and discussed. If you can't find your problem there, don't hesitate to drop me an E-Mail or open an issue yourself.
I am not responsible for any results produced with ViralClust nor for the conclusions you draw from it.

***

### Installation
In order to run ViralClust, I recommend creating a conda environment dedicated for NextFlow.
Of course, you can install NextFlow on your system how ever you like, but considering potential users not having sudo permissions, the conda way proofed to be simple.

* First install [conda](https://docs.conda.io/en/latest/) on your system: You'll find the latest installer [here](https://docs.conda.io/en/latest/miniconda.html).
* Next, make sure that conda is part of your `$PATH` variable, which is usually the case. For any issues, please refer to the respective [installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

<sub>**Warning**: Currently, ViralClust is running on Linux systems. Windows and MacOS support may follow at some point, but has never been tested so far.</sub>
* Create a conda environment and install NextFlow within this environment:

  <details><summary>Click here to see how:</summary>

  ```bash
  conda create -n nextflow -c bioconda nextflow
  conda activate nextflow
  ```
  </details>

  <sub>**Alternative**: You can also use the `environment.yml` provided, after cloning the repo:</sub>

  <details><summary>Click here to see how:</summary>

  ```bash
  conda env create -f environment.yml
  ```
   </details>

* Clone the github repository for the latest version of ViralClust, or download the latest stable release version here.

  <details><summary>Click here to see how:</summary>

  ```bash
  `git clone https://github.com/klamkiew/viralClust.git && cd viralClust`
  ```
   </details>

* Done!

***

### Quickstart

You may ask yourself, how you can run ViralClust yourself now.
Well, first of all, you do not have to worry about any dependencies, since NextFlow will take care of this via individual conda environments for each step. You just have to make sure to have a stable internet connection, when you run ViralClust for the very first time.
If not done yet, now is the time to activate your `conda` environment:
`conda activate nextflow`
And we're ready to go!

`nextflow run viralclust.nf --fasta "data/test_input.fasta"`


This might take a little bit of time, since all individual `conda` environments for each step of ViralClust is created.
In the mean time, let's talk about parameters and options.

### Parameters & Options

Let us briefly go over the most important parameters and options. There is a more detailed overview of all possible flags, parameters and additional stuff you can
do in the help of message of the pipeline - and at the end of this file.

***

### Graphical Workflow

![Workflow graph](/pic/workflow.png)

### Help Message

This paragraph is simply the help message of ViralClust.

<details><summary>Expand here</summary>

```
____________________________________________________________________________________________

Welcome to ViralClust - your pipeline to cluster viral genome sequences once and for all!
____________________________________________________________________________________________

Usage example:
nextflow run viralclust.nf --update_ncbi

or

nextflow run viralclust.nf --fasta "genomes.fasta"

or both

nextflow run viralclust.nf --update_ncbi --fasta "genomes.fasta"

____________________________________________________________________________________________

Mandatory Input:
--fasta PATH                      Path to a multiple fasta sequence file, storing all genomes that shall be clustered.
                                  Usually, this parameter has to be set, unless the parameter --ncbi_update has been set.

Optional Input:
--goi PATH                        Path to a (multiple) fasta sequence file with genomes that have to end
                                  up in the final set of representative genomes, e.g. strains of your lab that are
                                  of special interest. This parameter is optional.
____________________________________________________________________________________________

Options:
--eval                            After clustering, calculate basic statistics of clustering results. For each
                                  tool, the minimum, maximum, average and median cluster sizes are calculated,
                                  as well as the average distance of two representative genomes.

--ncbi                            Additionally to the evaluation performed by --eval, NCBI metainformation
                                  is included for all genomes of the input set. Therefore, the identifier of fasta records are
                                  scanned for GenBank accession IDs, which are then used to retrieve information about the taxonomy,
                                  accession date and accession country of a sequence. Implicitly calls --eval.
                                  Attention: If no database is available at data, setting this flag
                                  implicitly sets --ncbi_update.

--ncbi_update                     Downloads all current GenBank entries from the NCBI FTP server and processes the data to
                                  the databank stored at data.

Cluster options:
--cdhit_params                    Additional parameters for CD-HIT-EST cluster analysis. [default -c 0.9]
                                  For more information and options, we refer to the CD-HIT manual.

--hdbscan_params                  Additional parameters for HDBscan cluster analysis. [default ]
                                  For more information and options, please use
                                  nextflow run viralclust.nf --hdbscan_help or python3 bin/hdbscan_virus.py -h.

--sumaclust_params                Additional parameters for sumaclust cluster analysis. [default ]
                                  For more information and options, we refer to the sumaclust manual.

--vclust_params                   Additional parameters for vsearch cluster analysis. [default --id 0.9]
                                  For more information and options, we refer to the vsearch manual.

--mmseqs_params                   Additional parameters for MMSeqs2 cluster analysis. [default ]
                                  For more information and options, we refer to the MMSeqs2 manual.

Computing options:
--cores INT                       max cores per process for local use [default 1]
--max_cores INT                   max cores used on the machine for local use [default ALL]
--memory INT                      max memory in GB for local use [default 16.GB]
--output PATH                     name of the result folder [default viralclust_results]
--permanentCacheDir PATH          location for auto-download data like databases [default data]
--condaCacheDir PATH              location for storing the conda environments [default conda]
--workdir PATH                    working directory for all intermediate results [default /tmp/nextflow-work-$USER]

Nextflow options:
-with-report rep.html             cpu / ram usage (may cause errors)
-with-dag chart.html              generates a flowchart for the process tree
-with-timeline time.html          timeline (may cause errors)
____________________________________________________________________________________________

```
</details>