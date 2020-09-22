# Bioinformatic methods for Giagacomo et al 2020, Comparing DNA Extraction and 16s Amplification Methods for Plant-Associated Bacterial Communities 

This repo contains all the bioinformatic scripts necessary to run the analyses in this paper (currently in submission). 

## Folders
* **0_Troubleshooting** - Exactly what it says: minor scripts to do troubleshooting or otherwise test things out over the course of the analysis. Some of this is preliminary data, and all is only tangentally related to the main paper analyses
* **Figures** - SVG versions of the figures from the published manuscript
* **MakeBlockingOligos** - The support scripts for creating the maize-specific blocking oligos used in the paper
* **TestExtraction** - All scripts related to Experiment I in the paper, testing different DNA extraction methods
* **TestPrimers** - All scripts related to Experiment II in the paper, testing different methods of blocking plastid amplification
* **(root)** - The root folder contains the final analyses used in the paper. These are refinements of ones in the TestExtraction and TestPrimers folders, generally to put them into the final format for publication.

## How to run
All pipelines in this analysis are coordinated by a series of Bash scripts. All Bash scripts (and most support Python & R scripts) are prefaced with a number to indicate their position in the pipeline. (So, **1c_PlotPrimersInSamples.py** will be run before **1g_GraphDeblurReadRetention.py**, and there are presumably several steps in-between in the Bash script.)

To recreate our analyses, first run **0_CreateCondaEnvironments.sh** (in the root directory) to create a Conda environment with the necessary software for all analyses. (You must have some version of Conda installed for this to work; we use [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

If you want to run the individual analyses from **TestExtraction** and **TestPrimer**, go into those directories and start with the **0_Analyze16sExtractionMethods.sh** or **0_Analyze16sAmpMethods.sh** scripts, respectively, which are the master scripts that call all others. The necessary variables (data locations, read depth, etc.) are all defined inside these scripts and passed as arguments to downstream scripts. (Note that rerunning from scratch will require downloading the [raw data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA646931) from NCBI.)

To regenerate just the figures & tables from the manuscript, run **1_MakePublicationGraphics.sh** in the root directory, since all support files necessary to run this have been included in the repo.

## Software
The following software packages were used in this analysis (see [the paper](https://www.biorxiv.org/content/10.1101/2020.07.23.217901v1) for full citations).  All analyses were done on an desktop workstation with a 4-core Intel Xenon W-2123 processor and 64 GB of RAM running Linux Mint 18.3.

**R packages**
* ape v5.3 (Paradis and Schliep 2019)
* argparse v2.0.1 (Davis 2019)
* DESeq2 v1.24.0 (Love, Huber, and Anders 2014)
* dplyr v0.8.3 (Wickham et al. 2019)
* ggplot2 v3.2.1(Wickham 2016)
* gridExtra v2.3 (Auguie 2017)
* igraph v1.2.5 (Csardi, Nepusz, and Others 2006)
* phyloseq v1.28.0 (McMurdie and Holmes 2013)
* rbiom v1.0.0 (Smith 2019)
* tidyr v1.0.0 (Wickham and Henry 2019)
* vegan v2.5.5 (Oksanen et al. 2019)
  
**Python packages**
* argparse v1.1
* biopython v1.74 (Cock et al. 2009)
* matplotlib 3.1.0 (Hunter 2007)
* primer3-py v0.6.0 (“primer3-Py” n.d.)

**Command-line tools**
* biom v2.1.7 (McDonald et al. 2012)
* blast+ v2.2.31 (Camacho et al. 2009)
* Conda v4.8.2, Clustal Omega v1.2.1 (Sievers et al. 2011)
* cutadapt v2.4 (Martin 2011)
* libprimer3 v2.5.0 (“Primer3” n.d.)
* QIIME2 v2019-7 (Bolyen et al. 2019)
* vsearch v2.7.0 (Rognes et al. 2016)
