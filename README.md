# Calculate gene percentage  [![DOI](https://zenodo.org/badge/87009573.svg)](https://zenodo.org/badge/latestdoi/87009573)
A perl script that calculates estimated gene percentage in genome assembly given a FASTA file containing list of genes.

## Prerequistes
* [Perl5](https://dev.perl.org/perl5/)
* [BioPerl](http://bioperl.org) - This script was tested using BioPerl-1.6.1

## Run
To calculate estimated gene coverage:
1. Create a blast database of the output assembly
2. Align genes sequences against the blast database 

This could be done using the following [blast](https://blast.ncbi.nlm.nih.gov/) commands:
```
makeblastdb -in <PATH/TO/Assembly.fa> -out databaseBLAST -dbtype nucl -parse_seqids
blastn -query <PATH/TO/gene_file.fa> -out output.blast.txt -db databaseBLAST -num_threads <N>
```
3. Run gene_coverage.pl
```
perl gene_coverage.pl output.blast.txt
```
