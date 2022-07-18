#!/usr/bin/env python3

import os
import sys
from Replidec.Replidec_multi import bayes_classifier_batch,bayes_classifier_contig,bayes_classifier_genomes
from argparse import RawTextHelpFormatter,ArgumentParser


current_work_dir = os.path.dirname(os.path.realpath(__file__))

parser = ArgumentParser(description="replication cycle detector", formatter_class=RawTextHelpFormatter)
parser.add_argument('--version', action='version', version='Replidec v0.2.1')

parser.add_argument("-p","--program", default='multiSeqEachAsOne', required=True,
                    choices=['multiSeqAsOne','batch','multiSeqEachAsOne'],
                    help="multiSeqAsOne mode: input is a plain text file contain two coloumn (seprator must be **tab**)\n"
                        "   first column: sample name;\n"
                        "   second column: path of the protein file from one virues;\n"
                        "multiSeqEachAsOne mode: input is a sequence file and treat each seqence as from one virus;\n"
                        "batch mode: input is a plain text file contain two coloumn (seprator must be **tab**);\n"
                        "   first column: sample name;\n"
                        "   second column: path of the protein file from one virues;\n")

parser.add_argument("-i", "--input_file",
                    help="input file. Can be a sequence file or index file\n"
                        "multiSeqAsOne mode: input is sequence file\n"
                        "multiSeqEachAsOne mode: input is a sequence file\n"
                        "batch mode: input is a plain text file contain two coloumn (seprator must be **tab**)\n"
                        "   first column: sample name;\n"
                        "   second column: path of the protein file from one virues;\n")
parser.add_argument("-w", "--wd",default="Replidec", dest="workdir", help="work dir path")

parser.add_argument("-s", "--summary",default="BC_predict.summary", help="name of summary file")

parser.add_argument("-t", "--threads",default=10, type=int, help="number of parallel threads")

parser.add_argument("-c", "--hc",default=1e-5, type=float, dest="hmmer_creteria", help="Creteria to filter hmmer result")

parser.add_argument("-H", "--hp",default="--noali --cpu 3", dest="hmmer_parameter", help="Parameter used for hmmer")

parser.add_argument("-m", "--mc",default=1e-5, type=float, dest="mmseqs_creteria", help="Creteria to filter mmseqs2 result")

parser.add_argument("-M", "--mp", dest="mmseqs_parameter",
                    default="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3",
                    help="Parameter used for mmseqs")

parser.add_argument("-b", "--bc",default=1e-5, type=float, dest="blastp_creteria", help="Creteria to filter blast result")

parser.add_argument("-B", "--bp",default="-num_threads 3", dest="blastp_parameter", help="Parameter used for blastp")

args = parser.parse_args()


def main():
    print("Using %s"%args.program)
    if args.program == "multiSeqAsOne":
        bayes_classifier_genomes(args.input_file, args.workdir, summaryfile=args.summary, threads=args.threads,
                                hmm_creteria=args.hmmer_creteria, mmseqs_creteria=args.mmseqs_creteria, blastp_creteria=args.blastp_creteria,
                                hmmer_para=args.hmmer_parameter, mmseqs_para=args.mmseqs_parameter, blastp_para=args.blastp_parameter)
    elif args.program == "multiSeqEachAsOne":
        bayes_classifier_contig(args.input_file, args.workdir, summaryfile=args.summary, threads=args.threads,
                               hmm_creteria=args.hmmer_creteria, mmseqs_creteria=args.mmseqs_creteria, blastp_creteria=args.blastp_creteria,
                               hmmer_para=args.hmmer_parameter, mmseqs_para=args.mmseqs_parameter, blastp_para=args.blastp_parameter)
    elif args.program == "batch":
        bayes_classifier_batch(args.input_file, args.workdir, summaryfile=args.summary, threads=args.threads,
                              hmm_creteria=args.hmmer_creteria, mmseqs_creteria=args.mmseqs_creteria, blastp_creteria=args.blastp_creteria,
                              hmmer_para=args.hmmer_parameter, mmseqs_para=args.mmseqs_parameter, blastp_para=args.blastp_parameter)
    #elif args.program == "test_multiSeqAsOne":
    #    bayes_classifier_genomes("%s/example/genome_test.index"%current_work_dir, "./test_multiSeqAsOne",
    #                            summaryfile="BC_predict.summary",threads=10)
    #elif args.program == "test_multiSeqEachAsOne":
    #    bayes_classifier_contig("%s/example/test.contig.small.fa"%current_work_dir, "./test_multiSeqEachAsOne",
    #                            summaryfile="BC_predict.summary",threads=10)
    #elif args.program == "test_batch":
    #    bayes_classifier_batch("%s/example/example.list"%current_work_dir, "./test_batch",
    #                            summaryfile="BC_predict.summary")
    else:
        print("Please check the vaild program (multiSeqAsOne|multiSeqEachAsOne|batch)")

#### main ###
if __name__ == "__main__":
    main()
