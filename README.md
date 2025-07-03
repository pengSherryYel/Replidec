# Replidec: Replication Cycle Decipher for Phages

[![PyPI](https://img.shields.io/pypi/v/Replidec.svg)](https://pypi.python.org/pypi/Replidec)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/replidec/badges/version.svg)](https://anaconda.org/bioconda/replidec)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/replidec/badges/downloads.svg)](https://anaconda.org/bioconda/replidec)

## Aim

Use bayes classifier combine with homology search to predict virus replication cycle

## Install

### Method 1: using Conda (Recommend using bioconda with latest version)

```bash
conda create -n replidec
conda activate replidec
conda install -c conda-forge -c bioconda replidec
or
conda install -c denglab -c conda-forge -c bioconda replidec
```

### Method 2: using Docker

```bash
docker pull quay.io/biocontainers/replidec:0.3.5--pyhdfd78af_0
docker run quay.io/biocontainers/replidec:0.3.5--pyhdfd78af_0 Replidec -h
## Example
docker run -v /your/host/data:/data/ quay.io/biocontainers/replidec:0.3.5--pyhdfd78af_0 Replidec -i data/your_inputfile -p
choose_mode_based_on_your_input_type -w data
```

### Method 3: using pip

If you install using pip, please make sure that `mmseqs`, `hmmsearch` and `blastp` is set to $PATH, these software can equal or higher than version list below

- MMseqs2 Version: 13.45111

- HMMER 3.3.2 (Nov 2020)

- Protein-Protein BLAST 2.5.0+

```bash
pip3 install Replidec
```

## Usage: Overview

```
Replidec, Replication cycle prediction tool for prokaryotic viruses

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -p , --program        { multi_fasta | genome_table | protein_table }
                        
                        multi_fasta mode:
                        input is a fasta file and treat each sequence as one virus
                        
                        genome_table mode:
                        input is a tab separated file with two columns
                        ___1st column: sample name
                        ___2nd column: path to the genome sequence file of the virus
                        
                        protein_table mode:
                        input is a tab separated file with two columns
                        ___1st column: sample name
                        ___2nd column: path to the protein file of the virus
                        
  -i , --input_file     The input file, which can be a sequence file or an index table
  -w , --work_dir       Directory to store intermediate and final results (default = ./Replidec_results)
  -n , --file_name      Name of final summary file (default = prediction_summary.tsv)
  -t , --threads        Number of parallel threads (default = 10)
  -e , --hmmer_Eval     E-value threshold to filter hmmer result (default = 1e-5)
  -E , --hmmer_parameters 
                        Parameters used for hmmer (default = --noali --cpu 3)
  -m , --mmseq_Eval     E-value threshold to filter mmseqs2 result (default = 1e-5)
  -M , --mmseq_parameters 
                        Parameter used for mmseqs
                        (default = -s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3)
  -b , --blastp_Eval    E-value threshold to filter blast result (default =1e-5)
  -B , --blastp_parameter 
                        Parameters used for blastp (default = -num_threads 3)
  -d, --db_redownload   Remove and re-download database
```

## Usage: Download database (-d)

Database used in Replidec will be download automatically. 

Location: will be download at the where Replidec installed

If you want to redownload the database, `-d` parameter can be used. The older database will be mv to "discarded_db" in the workdir(-w); This dir can be removed manually by user.


## Usage: Input (-i) and Propgram (-p)

**Input file is different base on different program**

Replidec cantain **3** different program:

1. 'multi_fasta'
2. 'genome_table'
3. 'protein_table',

### multi_fasta mode:
* input is a **fasta** file and treat each sequence as one virus.
  * Example: <your_path>/viral_contigs.fasta
    
    ```
    >contig_1
    TATCGATCGATCGATCGATCGATCGTACGTACGTACGTACG...
    >contig_2
    CATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
    ...
    ```

### genome_table mode:
* input is a **tab** separated file with two columns.
    
    * 1st column: sample name
    * 2nd column: path to the genome sequence file of the virus
    * Example: <your_path>/example_genomes.tsv
 
    ```
    contig_1    your/file/path/contig_1.fasta
    contig_2    your/file/path/contig_2.fasta
    contig_3    your/file/path/contig_3.fasta
    ...
    ```
  
### protein_table mode:
* input is a **tab** separated file with two columns

    * 1st column: sample name
    * 2nd column: path to the protein file of the virus
    * Example: <your_path>/example_proteins.tsv

    ```
    contig_1_prot	your/file/path/contig_1.fasta
    contig_2_prot	your/file/path/contig_2.fasta
    contig_3_prot   your/file/path/contig_3.fasta
    ...
    ```

## Usage: Output (-w and -n)
The output directory can be assigned with `-w , --work_dir ` where the intermidiate files and the final prediction results will be stored.
The name of the final summary file can be assigned with `-n , --file_name` argument.

At the end of the analysis, the output directory would contain the following:
* BC_Inno: This directory contains the result file for dectect Innovirues
* BC_mmseqs: This directory contains the result file for mapping result to our custom database
* BC_pfam: This directory contains the result file for dectect the Integrase and Excisionase
* BC_prodigal: This directory contains the result file for CDS prediction from genome or contig sequence. (if {-p protein_table} is used, this directory will not be created)
* prediction_summary.tsv: This file is the summary file of the predict result. It contain multiple coloumns.
    * sample_name: identifier. Can be sequence id or first coloumn the plain text input file. 

    * integrase_number: the number of genes mapped to integrase meet the creteria(set by -c).

    * excisionase_number: the number of genes mapped to excisionase meet the creteria(set by -c).

    * pfam_label: if contain integrase or excisionase, label will be "Temperate". otherwise "Virulent".

    * bc_temperate: conditional probability of temperate|genes. 

    * bc_virulent: conditional probability of virulent|genes. 

    * bc_label: if bc_temperate greater than bc_virulent, label will be "Temperate". otherwise "Virulent".

    * final_label: if pfam_label and bc_label both is Temperate, then label will be "Temperate"; if Innovirues marker gene exist, then label will be "Chronic"; otherwise "Virulent".

    * match_gene_number:  the number of genes mapped to our custom databse.

    * path: path of input faa file


## Example (Data in test folder, please navigate to test folder first)
```
cd test

## Conda
## test passed - genome_table
replidec -p genome_table -i example/genome_test.small.index -w opt_folder_genome_table

## test passed - multi_fasta
replidec -p multi_fasta -i example/test.contig.small.fa -w opt_folder_multi_fasta

## test passed - protein_table
replidec -p protein_table -i example/example.small.list -w opt_folder_protein_table


## Docker
docker run -v /Your_path_clone_replidec/Replidec/test:/data/ quay.io/biocontainers/replidec:0.3.5--pyhdfd78af_0 Replidec -p multi_fasta -i /data/example/test.contig.small.new.fa -w /data/opt_folder_docker_multi_fasta
```


## Issues
### Database can not be downloaded automatically 
If the dataset cannot be automatically downloaded from Zenodo due to regional access restrictions, you may manually add it instead. The same database has also been uploaded to OSF as an alternative source.

1. **Locate your Replidec installation path**  
After installing Replidec via Conda or Docker, locate the installed directory. Typically, it can be found at:
`your_conda_path/envs/env_name/lib/python*/site-packages/Replidec`


2. **Navigate to the Replidec folder**  
Use the terminal to move into the directory:  
`cd your_conda_path/envs/env_name/lib/python*/site-packages/Replidec`

3. **Download the database manually from OSF (Project name: Replidec)**  
Access the alternative download link here:
👉 https://osf.io/thpkb/files/osfstorage

4. **Extract the database**  
After downloading, extract the contents of the archive into the Replidec directory and a folder named "db" will be created:
`tar -zxvf db_v0.3.2.tar.gz`  
✅ Note: Make sure the extracted folder can be found in this path `your_conda_path/envs/env_name/lib/python*/site-packages/Replidec/db`.

For now, everything is fixed, enjoy play with replidec.



