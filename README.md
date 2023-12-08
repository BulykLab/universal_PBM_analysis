# universal_PBM_analysis pipeline for ERISTWO

Pipeline originally written by Luis Barrera, Katy Weinand, _et al_.

Adapted by Raehoon Jeong 12/08/2023

## 1. Description
This repository contains code to run the universal PBM analysis pipeline. It begins with .gpr files and outputs position weight matrices (PWMs), frequency matrices, motif logos, and an html report that aggregates the result, among other files. For a thorough description of the steps in the pipeline, also refer to this [overview](https://github.com/BulykLab/universal_PBM_analysis/blob/main/Overview.md). Currently, this directory is located in `/data/bulyk/pipelines/universal_PBM_analysis` in ERIS.


## 2. Setting up a virtual environment for the PBM analysis pipeline

### Step 0 Logging into ERISTWO
For further information about the ERISTWO computing cluster, refer to the information page (https://rc.partners.org/kb/article/1315). The code below takes you to an interactive node in ERISTWO. You should only run jobs in the interactive node or submit them to the compute node.

```
## 0-1. Log in to ERISTWO
ssh yourid@eristwo.partners.org

## 0-2. Start an interactive job
bsub -Is -q interactive -R 'rusage[mem=16000]' /bin/bash'
```

### Step 1 Creating a virtual environment
In order to make sure that all the prerequisite packages are installed with the correct version of python (3.7) and R (3.6), setting up a conda virtual environment is useful.

```
## 1-1. Create a virtual environment (I named it pbmenv)
conda create -n pbmenv python=3.7 R=3.6 -y

## 1-2. Activate the PBM analysis virtual environment (pbmenv)
conda activate pbmenv
```

### Step 2 Installing necessary packages
Here is a list of specific packages that need to be installed: seqlogo, minpack, and weblogo. However, if there are other packages that are missing, use the `conda install [package]` command to install that package.

```
conda install -c conda-forge numpy -y
conda install -c bioconda bioconductor-seqlogo -y
conda install -c conda-forge r-minpack.lm -y
conda install -c conda-forge weblogo -y
```

## 3. Test run
Try running this test to make sure that the virtual environment is set up properly and that the necessary packages are installed. The user must create a directory within the `testfile` directory and copy the gpr files to run. Also, 

```
## 0-1. Log in to ERISTWO
ssh yourid@eristwo.partners.org

## 0-2. Start an interactive job
bsub -Is -q interactive -R 'rusage[mem=16000]' /bin/bash'

## 0-3. Move to PBM analysis pipeline directory
cd /data/bulyk/pipelines/universal_PBM_analysis

## 1. Activate the PBM analysis virtual environment (pbmenv)
## YOU MUST RUN THIS BEFORE ANY JOB SUBMISSION!
conda activate pbmenv

## 2. Create a directory [dirname] with test data (_i.e._ gpr files)
### This directory is where the input gpr files and the output files 
cd testfile 
mkdir [dirname]
cp test_data/*gpr [dirname]/

## 3. Copy and modify the directory name to [dirname] in the test job submission code
cp test_pbm_pipeline.lsf [filename].lsf
## Open the file and change line 27 to 'DATA=[dirname]'
vim [filename].lsf

## 4. Submit the job
bsub < [filename].lsf

## 5. Check the result in [dirname]. Open report.html to see if the result looks as expected.
```
