#!/usr/bin/bash
## test passed - multiSeqAsOne
Replidec -p multiSeqAsOne -i example/genome_test.small.index -w multiSeqAsOne
## test passed - multiSeqEachAsOne
Replidec -p multiSeqEachAsOne -i example/test.contig.small.fa -w multiSeqEachAsOne
## test passed - batch
Replidec -p batch -i example/example.small.list -w batch
