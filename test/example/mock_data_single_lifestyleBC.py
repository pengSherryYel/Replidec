# coding: utf-8
# %load mock_data_single_lifestyleBC.py
import sys
#sys.path.append("/home/viro/xue.peng/script/bayes_classifier/")
sys.path.append("/home/viro/xue.peng/script/bayes_classifier_v2/")
#from lifestyleBC import bayes_classifier_single,bayes_classifier_batch,bayes_classifier_contig
from Replidec import bayes_classifier_single,bayes_classifier_batch,bayes_classifier_contig
#bayes_classifier_batch("./sample1.acc.path","./batch_sample1.single","BC_predict.sample1.single.summary")
#bayes_classifier_batch("./sample2.acc.path","./batch_sample2.single","BC_predict.sample2.single.summary")
#bayes_classifier_batch("./sample3.acc.path","./batch_sample3.single","BC_predict.sample3.single.summary")

#infa=sys.argv[1]
#name=sys.argv[2]
infa="./test.contig.fa"
name="test.contig"
bayes_classifier_contig(infa,name,"%s.summary"%name)
