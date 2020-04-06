This directory contains relative abundance data from a species abundance profile table for 250 human feces samples fetched from https://db.cngb.org/microbiome/genecatalog/genecatalog/?gene_name=Human%20Gut%20(11M) </br>
Species from the table were searched on NCBI's Nucleotide database with the schema: </br>
result=entrez_search(db="nuccore", term=paste(trueID[x], "whole genome", "AND refseq[filter]", "AND 500000:1000000000[Sequence Length]", "NOT \"sequencing project\"",  sep=" ")) </br>
Any species whose genomes that could not be fetched were removed from the table and relative abundance values were scaled to equal 1
