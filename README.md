# MR_BETA
Scripts associated with the tool MR BETA: Metagenomic Reassignment based on Bias Estimation of Taxonomic Abundance
</br> </br>
## Tutorial </br>
The following steps have been designed to run on Pegasus, a high performace computing environment available at The George Washington University. Similar steps can likely be followed on most high performance computing environments that use bash. </br> </br>
### Step 1: </br>
Install miniconda in your local directory: https://docs.conda.io/en/latest/miniconda.html </br>
Miniconda will allow you to seemlessly install other packages you will need. It will also help manage dependencies. </br> </br>
### Step 2: </br>
Make sure that `bioconda` is included in your list of channels. To do so, open the file `.condarc` in a text editor of your choice (ie. nano, vim) and add `- bioconda` in your list of channels. Your .condarc should look something like this:
```channels:
  - bioconda
  - conda-forge
  - defaults
```
### Step 2: <br>
Install git by running `conda install git` </br> </br>

### Step 3: </br>
Clone BEAR into your your local directory by running `git clone https://github.com/sej917/BEAR.git` </br>
BEAR (Better emulation for artificial reads) is a tool developed by Stephen Johnson, Brett Trost, Dr. Jeffrey R. Long, and Dr. Anthony Kusalik of the University of Saskatchewan, Department of Computer Science. </br> </br>

### Step 4: </br>
Install Kraken 2 in a coda environment called *kraken2* by running `conda create -n kraken2 kraken2` </br> </br>

### Step 5:
