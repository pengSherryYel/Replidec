# Replidec: Replication Cycle Decipher for Phages

[![PyPI](https://img.shields.io/pypi/v/Replidec.svg)](https://pypi.python.org/pypi/Replidec)
[![Anaconda-Server Badge](https://anaconda.org/denglab/replidec/badges/version.svg)](https://anaconda.org/denglab/replidec)
[![Anaconda-Server Badge](https://anaconda.org/denglab/replidec/badges/license.svg)](https://anaconda.org/denglab/replidec)
[![Anaconda-Server Badge](https://anaconda.org/denglab/replidec/badges/downloads.svg)](https://anaconda.org/denglab/replidec)

## Aim

Use bayes classifier combine with homology search to predict virus replication cycle

## Install

### Method 1: using Conda

```bash
conda create -n replidec
conda activate replidec
conda install -c conda-forge -c bioconda replidec
or
conda install -c denglab -c conda-forge -c bioconda replidec
```

### Method 2: using Docker

```bash
docker pull denglab/replidec
```

If you want to use `Replidec` on an HPC, singularity is recommended. You can create a singularity image using following command,

```bash
singularity pull replidec.sif docker://denglab/replidec
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
Replidec [-h] [--version] -p {multiSeqAsOne,batch,multiSeqEachAsOne}
         [-i INPUT_FILE] [-w WORKDIR] [-s SUMMARY] [-t THREADS] [-c HMMER_CRETERIA] [-H HMMER_PARAMETER] [-m MMSEQS_CRETERIA]
         [-M MMSEQS_PARAMETER] [-b BLASTP_CRETERIA] [-B BLASTP_PARAMETER] [-d] [-D]
```
## Usage: database (-d & -D)

Database used in Replidec will be download automatically. 

Location: will be download at the where Replidec installed

If you want to redownload the database, `-d` parameter can be used. The older database will be mv to "discarded_db" in the workdir(-w); This dir can be removed manually by user.

We provide 2 database for user choose. Use -D to choose database. all|prokaryote

1. all: This database contain all proteins from the prokaryotic, eukaryotic
virus and prophage.
2. prokaryote: This database contain proteins from the prokaryotic virus and
prophage. eukaryotic virus gene was removed from this database

## Usage: Input(-i) and Propgram(-p)

**Input file is different base on different program**

Replidec cantain **3** different program:

1. 'multiSeqAsOne'
2. 'batch'
3. 'multiSeqEachAsOne',

### multiSeqAsOne
* multiSeqAsOne mode: input is a plain text file contain two coloumn (seprator must be **tab**)

    * first column: sample name; this will be used as identfier in the output summary file 
    
    * second column: path of the **genome or contig** file from one virues (Each file can contain multi seq)

    * Example: test/example/genome_test.small.index

    ```
    seq1    example/genome_test/genome.test.fnaaa
    seq2    example/genome_test/genome.test.fnaab
    seq3    example/genome_test/genome.test.fnaac
    ```

### multiSeqEachAsOne
* multiSeqEachAsOne mode: input is a **sequence** file and treat *each* seqence as from one virus and give each sequence a predict result;
    
    * This mode will treat each sequence independently

    * Example: test/example/test.contig.small.fa

### batch
* batch mode: input is a plain text file contain two coloumn (seprator must be **tab**);

    * first column: sample name;

    * second column: path of the **protein** file from one virues;

    * Example: test/example/example.small.list

    ```
    simulate_art_sample1.10 example/simulate_art_sample1.10.faa
    simulate_art_sample1.11 example/simulate_art_sample1.11.faa
    simulate_art_sample1.12 example/simulate_art_sample1.12.faa
    ```

## Usage: Output(-w and -s)
The output dirname can use `-w` to set and the name of summary file can use `-s` to set.
Under output dir serveral dir and a summary file will be generated
* BC_Inno: This dir contain the result file for dectect Innovirues
* BC_mmseqs: This dir contain the result file for mapping result to our custom database
* BC_pfam: This dir contain the result file for dectect the Integrase and Excisionase
* BC_prodigal: This dir contain the result file for CDS prediction from genome or contig sequence. (-p batch will not generate this dir)
* BC_predict.summary: This file is the summary file of the predict result. It contain multiple coloumns.
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

## Example
```
## test passed - multiSeqAsOne
Replidec -p multiSeqAsOne -i example/genome_test.small.index -w multiSeqAsOne

## test passed - multiSeqEachAsOne
Replidec -p multiSeqEachAsOne -i example/test.contig.small.fa -w multiSeqEachAsOne

## test passed - batch
Replidec -p batch -i example/example.small.list -w batch
```

