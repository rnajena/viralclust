# WORK IN PROGRESS
There are some bugs known to me, some aren't. Before getting desperate, please check out the Issues that are already opened and discussed. If you can't find your problem there, don't hesitate to drop me an E-Mail or open an issue yourself.
I am not responsible for any results produced with ViralClust nor for the conclusions you draw from it.

***
## ViralClust - Find representative viruses for your dataset
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![Python3.6](https://img.shields.io/badge/Language-Python_3.6-teal.svg)[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/klamkiewicz?label=%40klamkiewicz&style=social)](https://twitter.com/klamkiewicz)

***

### Quickstart


Clone the repository.

`git clone https://github.com/klamkiew/viralclust.git`

Create a conda environment for nextflow.

`conda create --name nextflow`

`conda activate nextflow`

and actually install nextflow.

`conda install -c bioconda nextflow`


Run this command to perform the (currently) full analyses: apply all cluster algorithm, build a MSA and a tree.

`nextflow run viralclust.nf --fasta "<INPUT FASTA>" --max_cores "<CPUs TO USE>" --workdir "<TMP-DIRECTORY>" --output "<OUTPUT DIRECTORY>" --tree`

---

![Workflow graph](/pic/flowchart.png)