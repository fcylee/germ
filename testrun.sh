#!/bin/bash

# Test run for standard sequences
Rscript germs.R -f test_data/test.fasta -o output.tsv.gz -w 122 -s 122 -t transcripts.txt -p plots

# Test run for sequences non-standard characters
Rscript germs.R -f test_data/test_nonstdchars.fasta -o output_nonstdchars.tsv.gz -w 122 -s 122 -t transcripts_nonstdchars.txt -p plots_nonstdchars

# Test run for standard sequences (no output filename provided)
Rscript germs.R -f test_data/test.fasta -w 122 -s 122 -t transcripts.txt

# Test run for pmd
Rscript germs.R -f test_data/test_pdm.fasta -d test_data/test_matrix.csv -w 122 -s 122 -t "2:29675391-29679862(-)" -p plots_pdm

# Test run for pmd comparing no matrix
Rscript germs.R -f test_data/test_pdm.fasta -w 122 -s 122 -t "2:29675391-29679862(-)" -p plots