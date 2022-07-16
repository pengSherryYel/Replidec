#!/usr/bin/bash
## test passed - test_multiSeqAsOne
#python ./Replidec_cmdline.py -p test_multiSeqAsOne
## test passsed - test_multiSeqEachAsOne
#python ./Replidec_cmdline.py -p test_multiSeqEachAsOne
## test passsed - test_batch
#python ./Replidec_cmdline.py -p test_batch
## test passed - multiSeqAsOne
#./Replidec_cmdline.py -p multiSeqAsOne -i example/genome_test.small.index
## test passed - multiSeqEachAsOne
#./Replidec_cmdline.py -p multiSeqEachAsOne -i example/test.contig.small.fa -w multiSeqEachAsOne
## test passed - batch
./Replidec_cmdline.py -p batch -i example/example.small.list -w batch
