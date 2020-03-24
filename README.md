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
```
channels:
  - bioconda
  - conda-forge
  - defaults
```
### Step 3: <br>
Install git by running `conda install git` </br> </br>

### Step 4: </br>
Clone this repository by running `git clone https://github.com/tyden46/MR_BETA.git` </br> </br>

### Step 5: </br>
Clone BEAR into your your local directory by running `git clone https://github.com/sej917/BEAR.git` </br>
BEAR (Better emulation for artificial reads) is a tool developed by Stephen Johnson, Brett Trost, Dr. Jeffrey R. Long, and Dr. Anthony Kusalik of the University of Saskatchewan, Department of Computer Science. </br> </br>

### Step 6: </br>
Install Kraken 2 in a coda environment called *kraken2* by running `conda create -n kraken2 kraken2` </br> </br>

### Step 7: </br>
Install ncbi-acc-download by running `conda create -n ncbi-acc-download ncbi-acc-download` </br>

### Step 8: </br>
Activate your ncbi-acc-download conda environment by running `conda activate ncbi-acc-download` </br> </br>

### Step 9: </br>
Fetch your genomes of interest by running `ncbi-acc-download --format fasta --out myGenomes.fna "$(< refSeqIds.csv)"` </br>
Note: This tutorial analysis uses a .csv formatted list of refSeqIds that constitute a "core" microbiome of eight healthy patients. In your own analysis, feel free to swap `refSeqIds.csv` for your own .csv formatted list of RefSeq IDs.

### Step 10: </br>
BEAR requires python 2.7. Create a conda environment called python2.7 by running `conda create -n python2.7 python=2.7` </br>
Run c`conda activate python2.7` </br>
Run BEAR via the command `python BEAR_release/scripts/generate_reads2.py -r myGenomes.fna -a coreMedianAbundance.txt -o 10Million-2.fna -t 10000000 -l 200`
