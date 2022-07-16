#!/usr/bin/bash
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
