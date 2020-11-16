
## <samp>ViralClust</samp> - Find representative viruses for your dataset
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.8](https://img.shields.io/badge/Language-Python_3.8-darkred.svg)![NextFlow](https://img.shields.io/badge/Nextflow-20.07.01-blue.svg)![conda](https://img.shields.io/badge/Uses-conda-green.svg)

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)

***

### DISCLAIMER
This pipeline is work-in-progress.
There are some bugs known to me, some aren't. Before getting desperate, please check out the Issues that are already opened and discussed. If you can't find your problem there, don't hesitate to drop me an E-Mail or open an issue yourself.
I am not responsible for any results produced with <samp>ViralClust</samp> nor for the conclusions you draw from it.

***

### Overview: What is this about?
Have you ever been in the situation that you wanted to compare you're specific virus of interest with all other viruses from its genus? Or even family? For some taxonomic clades, there are many different genomes available, which can be used for comparative genomics.

However, more often than not viral genome datasets are redundant and thus introduce bias into your downstream analyses. Think about a consensus genome of *Flavivirus* with 2.000 Dengue virus genomes and 5 Zika virus genomes. You may start to see the problem here. To remove redundancy, clustering of the input sequences is a nice idea. However, given the scientific question at hand, it is hard to determine whether a cluster algorithm is appropiate.

Thus, <samp>ViralClust</samp> was developed. A Nextflow pipeline utilizing different cluster methods and implementations all at once on your data set. Combining this with meta information from the NCBI allows you to explore the resulting representative genomes for each tool and decide for the cluster that fit your question.

For example: clustering all available *Filoviridae* with <samp>cd-hit-est</samp> usually leads to a large cluster containing all *Zaire Ebola viruses*, which can be valueable, if you want to compare this species as a whole. If you are interested in subtle changes within the species, you may want to use another approach, which divides the "Zaire cluster" into smaller sub-cluster, which represent different outbreaks and epidemics.

***

### Installation
In order to run <samp>ViralClust</samp>, I recommend creating a conda environment dedicated for NextFlow.
Of course, you can install NextFlow on your system how ever you like, but considering potential users not having sudo permissions, the conda-way proofed to be simple.

* First install [conda](https://docs.conda.io/en/latest/) on your system: You'll find the latest installer [here](https://docs.conda.io/en/latest/miniconda.html).
* Next, make sure that conda is part of your <samp>$PATH</samp> variable, which is usually the case. For any issues, please refer to the respective [installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

<sub>**Warning**: Currently, <samp>ViralClust</samp> is running on Linux systems. Windows and MacOS support may follow at some point, but has never been tested so far.</sub>
* Create a conda environment and install NextFlow within this environment:

  <details><summary>Click here to see how:</summary>

  ```bash
  conda create -n nextflow -c bioconda nextflow
  conda activate nextflow
  ```
  </details>

  <sub>**Alternative**: You can also use the <samp>environment.yml</samp> provided, after cloning the repo:</sub>

  <details><summary>Click here to see how:</summary>

  ```bash
  conda env create -f environment.yml
  ```
   </details>

* Clone the github repository for the latest version of <samp>ViralClust</samp>, or download the latest stable release version [here](https://github.com/klamkiew/viralClust/releases).

  <details><summary>Click here to see how:</summary>

  ```bash
  `git clone https://github.com/klamkiew/ViralClust.git && cd ViralClust`
  ```
   </details>

* Done!

***

### Quickstart

You may ask yourself, how you can run <samp>ViralClust</samp> yourself now.
Well, first of all, you do not have to worry about any dependencies, since NextFlow will take care of this via individual conda environments for each step. You just have to make sure to have a stable internet connection, when you run <samp>ViralClust</samp> for the very first time.
If not done yet, now is the time to activate your <samp>conda</samp> environment:
`conda activate nextflow`
And we're ready to go!

`nextflow run viralclust.nf --fasta "data/test_input.fasta"`


This might take a little bit of time, since all individual <samp>conda</samp> environments for each step of <samp>ViralClust</samp> is created.
In the mean time, let's talk about parameters and options.

### Parameters & Options

Let us briefly go over the most important parameters and options. There is a more detailed overview of all possible flags, parameters and additional stuff you can
do in the help of message of the pipeline - and at the [end of this file](#help-message).

###### Input sequences: --fasta \<PATH>
<samp>--fasta \<PATH\></samp> is the main parameter you **have** to set. This will tell <samp>ViralClust</samp> where your genomes are located. <samp>\<PATH\></samp> refers to a multiple fasta sequence file, which stores all your genomes of interest.

###### Specific genomes of interest: --goi \<PATH>
<samp>--goi \<PATH\></samp> is similar to the <samp>--fasta</samp> parameter, but the sequences stored in this specfic fasta file are your **g**enomes **o**f **i**nterest, or shortly GOI. Using this parameter tells <samp>ViralClust</samp> to include all genomes present in <samp>goi.fasta</samp> in the final set of representative sequences. You have a secret in-house lab-strain that is not published yet? Put it in your <samp>goi.fasta</samp>.

###### Evaluate and rate cluster: --eval and --ncbi
<samp>--eval</samp> and <samp>--ncbi</samp> are two parameters, that do more for you than just clustering. Since <samp>ViralClust</samp> is running several clustering algorithms, it can be hard to decide which one produced the most appriopate results. Worry not, since <samp>--eval</samp> is here to help you. Additionally to the clustering results, you'll get a brief overview of the clusters, that arose from the different algorithms. With <samp>--ncbi</samp> enabled, <samp>ViralClust</samp> further scans your genome identifiers (the lines in your fasta file starting with <samp>\></samp>) for GenBank accession IDs and uses them to retrieve further information from the NCBI about the taxonomy of the sequence, as well as accession date and country. Note that using <samp>--ncbi</samp> implicitly also sets <samp>--eval</samp>.

###### Update the NCBI metainformation database: --update_ncbi
<samp>--update_ncbi</samp> is used whenever you need to update the database of <samp>ViralClust</samp>. As soon as you run the pipeline with <samp>--ncbi</samp> enabled for the first time, this is done automatically for you. Each viral GenBank entry currently available from the NCBI is processed and for each entry, <samp>ViralClust</samp> stores the accession ID, taxonomy, accession date and accession country for future uses.

###### Specify the output path: --output \<PATH>
<samp>--output \<PATH\></samp> specifies the output directory, where all results are stored. Per default, this is a folder called <samp>ViralClust_result</samp> which will be created in the directory that you are currently in.

###### Determine the numbers of cores used: --cores and --max_cores
<samp>--max_cores</samp> and <samp>--cores</samp> determine how many CPU cores are used at maximum and how many cores are used for one individual process at maximum, respectively. The default values cause <samp>ViralClust</samp> to use all available cores, but for each individual step in the pipeline, only 1 core is used.

There are many more parameters, especially directly connected to the behaviour of Nextflow, which are not explained here. The main things are covered, for the rest, I refer to the [clustering section](#cluster-tools) and the [complete help message](#help-message) of <samp>ViralClust</samp>.


***

### Cluster Tools

Since <samp>ViralClust</samp> is nothing without the great work of awesome minds, it is only fair to give credit, where credit is due. Currently, five different approaches are used, to cluster input genomes. <samp>[CD-HIT](http://www.bioinformatics.org/cd-hit/cd-hit-user-guide)</samp>, <samp>[Sumaclust](https://git.metabarcoding.org/obitools/sumaclust/wikis/home/)</samp> and <samp>[vsearch](https://github.com/torognes/vsearch)</samp> all implement the same algorithmic idea, but with minor, subtle changes in their respective heuristics. I further utilize the clustering module of <samp>[MMSeqs2](https://github.com/soedinglab/MMseqs2)</samp>. And, last but not least, <samp>ViralClust</samp> implements a <samp>k-mer</samp> based clustering method, which is realized with the help of <samp>[UMAP](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html)</samp> and <samp>[HDBSCAN](https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html)</samp>.

For all tools, the respective manual and/or github page is linked. Firstly, because I think, all of those are great tools, which you are implicitly using by using <samp>ViralClust</samp>. And second, because <samp>ViralClust</samp> offers the possibility to set all parameters of all tools; therefore, if you need something very specific, you can check out the respective documentations.

And, in case of using any of the results provided by <samp>ViralClust</samp> in a scientific publication, I would be grateful to be cited. In my eyes, it is only fair that you not only cite <samp>ViralClust</samp>, but also the clustering method you ultimately decided for, even if <samp>ViralClust</samp> was assisting you in the decision.

<details><summary>Click here for all citations</summary>

  * CD-HIT:
    * `Weizhong Li & Adam Godzik, "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences". Bioinformatics, (2006) 22:1658-9`
    * `Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, (2012), 28 (23): 3150-3152`

  * sumaclust:
    * `Mercier C, Boyer F, Bonin A, Coissac E (2013) SUMATRA and SUMACLUST: fast and exact comparison and clustering of sequences. Available: http://metabarcoding.org/sumatra.`

  * vsearch:
    * `Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584`

  * MMSeqs2:
    * `Steinegger, M., Söding, J. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets". Nat Biotechnol 35, 1026–1028 (2017)`

  * UMAP & HDBscan: 
    * `McInnes, L, Healy, J, "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", ArXiv e-prints 1802.03426, 2018`
    * `L. McInnes, J. Healy, S. Astels, "hdbscan: Hierarchical density based clustering" In: Journal of Open Source Software, The Open Journal, volume 2, number 11. 2017`
</details>

***

### Graphical Workflow

![Workflow graph](/pic/workflow.png)

### Help Message

This paragraph is simply the help message of <samp>ViralClust</samp>.

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
                                  You can use nextflow run viralclust.nf --cdhit_help
                                  For more information and options, we refer to the CD-HIT manual.

--hdbscan_params                  Additional parameters for HDBscan cluster analysis. [default -k 7]
                                  For more information and options, please use
                                  nextflow run viralclust.nf --hdbscan_help.

--sumaclust_params                Additional parameters for sumaclust cluster analysis. [default -t 0.9]
                                  You can use nextflow run viralclust.nf --sumaclust_help.
                                  For more information and options, we refer to the sumaclust manual.

--vclust_params                   Additional parameters for vsearch cluster analysis. [default --id 0.9]
                                  You can use nextflow run viralclust.nf --vclust_help
                                  For more information and options, we refer to the vsearch manual.

--mmseqs_params                   Additional parameters for MMSeqs2 cluster analysis. [default --min-seq-id 0.9]
                                  You can use nextflow run viralclust.nf --mmseqs_help
                                  For more information and options, we refer to the MMSeqs2 manual.

Computing options:
--cores INT                       max cores per process for local use [default 1]
--max_cores INT                   max cores used on the machine for local use [default 8]
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