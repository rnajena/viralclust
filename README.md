
## ViralClust - Find representative viruses for your dataset
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.8](https://img.shields.io/badge/Language-Python_3.8-darkred.svg)![NextFlow](https://img.shields.io/badge/Uses-NextFlow_DSL2-blue.svg)![conda](https://img.shields.io/badge/Uses-conda-green.svg)
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
  <details><summary><sub>Click here to see how:</sub></summary>

  ```bash
  conda create -n nextflow -c bioconda nextflow
  conda activate nextflow
  ```
  </details>

  <sub>**Alternative**: You can also use the `environment.yml` provided, after cloning the repo:</sub>
  <details><summary><sub>Click here to see how:</sub></summary>

  ```bash
  conda env create -f environment.yml
  ```
   </details>

* Clone the github repository for the latest version of ViralClust, or download the latest stable release version here.
  <details><summary><sub>Click here to see how:</sub></summary>

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


***

### Graphical Workflow

![Workflow graph](/pic/workflow.png)