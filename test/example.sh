#!/usr/bin/bash
## test passed - genome_table
replidec -p genome_table -i example/genome_test.small.index -w opt_folder_genome_table
## test passed - multi_fasta
replidec -p multi_fasta -i example/test.contig.small.fa -w opt_folder_multi_fasta
## test passed - protein_table
replidec -p protein_table -i example/example.small.list -w opt_folder_protein_table
