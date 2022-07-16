# Replidec: Replication Cycle Detector for Phages
## Aim

Use bayes classifier combine with homology search to predict virus replication cycle

### Important: 

mmseqs, hmmsearch, blastp must set to $PATH, these software can equal or higher than version list below

    MMseqs2 Version: 13.45111

    HMMER 3.3.2 (Nov 2020)

    Protein-Protein BLAST 2.5.0+

```
conda create -n replidec
conda activate replidec
conda install -c bioconda,anaconda "python>=3.8" mmseqs2 hmmer blast
```

## Install
```
pip3 install Replidec
or
pip3 install https://github.com/pengSherryYel/Replidec_v0.2.1/releases/download/v0.2.1/Replidec-0.2.1-py2.py3-none-any.whl
```

## Usage
```
Replidec [para]
```

Replidec cantain 6 different program:

'multiSeqAsOne', 'batch', 'multiSeqEachAsOne',

'test_multiSeqAsOne','test_batch','test_multiSeqEachAsOne' 

* multiSeqAsOne mode: input is a plain text file contain two coloumn (seprator must be **tab**)

   first column: sample name;

   second column: path of the protein file from one virues

* multiSeqEachAsOne mode: input is a sequence file and treat each seqence as from one virus;

* batch mode: input is a plain text file contain two coloumn (seprator must be **tab**);

   first column: sample name;

   second column: path of the protein file from one virues;

* test_* mode: test for each prpgram

## Example
```
## test passed - test_multiSeqAsOne
Replidec -p test_multiSeqAsOne

## test passsed - test_multiSeqEachAsOne
Replidec -p test_multiSeqEachAsOne

## test passsed - test_batch
Replidec -p test_batch

## test passed - multiSeqAsOne
Replidec -p multiSeqAsOne -i example/genome_test.small.index -w multiSeqAsOne

## test passed - multiSeqEachAsOne
Replidec -p multiSeqEachAsOne -i example/test.contig.small.fa -w multiSeqEachAsOne

## test passed - batch
Replidec -p batch -i example/example.small.list -w batch
```

