#!/usr/bin/env python3

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from Replidec.Replidec_multi import bayes_classifier_batch,bayes_classifier_contig,bayes_classifier_genomes
from argparse import RawTextHelpFormatter,ArgumentParser,ArgumentDefaultsHelpFormatter
from Replidec.Replidec import checkdb_and_download
from Replidec.utility import mkdirs


current_work_dir = os.path.dirname(os.path.realpath(__file__))

parser = ArgumentParser(description="Replidec, Replication cycle prediction tool for prokaryotic viruses", formatter_class=RawTextHelpFormatter)
parser.add_argument("-v", "--version", action='version', version='Replidec v0.3.5')

parser.add_argument("-p","--program", default='multi_fasta', required=True,
                    choices=['multi_fasta','genome_table','protein_table'],
                    metavar="",
                    help=
                        "{ multi_fasta | genome_table | protein_table }"
                        "\n\n"
                        
                        
                        "multi_fasta mode:\n"
                        "input is a fasta file and treat each sequence as one virus\n"
                        "\n"
                        "genome_table mode:\n"
                        "input is a tab separated file with two columns\n"
                        "___1st column: sample name\n"
                        "___2nd column: path to the genome sequence file of the virus\n"
                        "\n"
                        "protein_table mode:\n"
                        "input is a tab separated file with two columns\n"
                        "___1st column: sample name\n"
                        "___2nd column: path to the protein file of the virus\n"
                        "\n"

                    )

parser.add_argument("-i", "--input_file",  metavar="",
                    help="The input file, which can be a sequence file or an index table\n")

parser.add_argument("-w", "--work_dir", default="Replidec_results", metavar="", dest="working_directory", help="Directory to store intermediate and final results (default = ./Replidec_results)")

parser.add_argument("-n", "--file_name", default="prediction_summary.tsv", metavar="", dest="file_name", help="Name of final summary file (default = prediction_summary.tsv)")

parser.add_argument("-t", "--threads", default=10, type=int, metavar="", help="Number of parallel threads (default = 10)")

parser.add_argument("-e", "--hmmer_Eval", default=1e-5, type=float, metavar="", dest="hmmer_Evalue_threshold", help="E-value threshold to filter hmmer result (default = 1e-5)")

parser.add_argument("-E", "--hmmer_parameters", default="--noali --cpu 3", metavar="", dest="hmmer_parameters", help="Parameters used for hmmer (default = --noali --cpu 3)")

parser.add_argument("-m", "--mmseq_Eval", default=1e-5, type=float, metavar="", dest="mmseqs_Evalue_threshold", help="E-value threshold to filter mmseqs2 result (default = 1e-5)")

parser.add_argument("-M", "--mmseq_parameters", dest="mmseqs_parameters", metavar="",
                    default="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3",
                    help="Parameter used for mmseqs\n"
                    "(default = -s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3)"
                    )

parser.add_argument("-b", "--blastp_Eval", default=1e-5, type=float, metavar="", dest="blastp_Evalue_threshold", help="E-value threshold to filter blast result (default =1e-5)")

parser.add_argument("-B", "--blastp_parameter", default="-num_threads 3", metavar="", dest="blastp_parameters", help="Parameters used for blastp (default = -num_threads 3)")

parser.add_argument("-d", "--db_redownload", action='store_true', default=False, dest="db_redownload", help="Remove and re-download database")


args = parser.parse_args()


def main():
    if args.db_redownload:
        print("Check db")
        fileDir = os.path.dirname(os.path.abspath(__file__))
        checkdb_and_download(fileDir,redownload=True)

    print("Using %s"%args.program)
    if not args.input_file:
        raise ValueError("Please provide an input file with [-i/ --input_file] argument.")

    if args.program == "genome_table":
        bayes_classifier_genomes(args.input_file, args.working_directory,
                                summaryfile=args.file_name, threads=args.threads,
                                hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold, blastp_criteria=args.blastp_Evalue_threshold,
                                hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters, blastp_para=args.blastp_parameters)
    elif args.program == "multi_fasta":
        bayes_classifier_contig(args.input_file, args.working_directory,
                               summaryfile=args.file_name, threads=args.threads,
                               hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold, blastp_criteria=args.blastp_Evalue_threshold,
                               hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters, blastp_para=args.blastp_parameters)
    elif args.program == "protein_table":
        bayes_classifier_batch(args.input_file, args.working_directory,
                               summaryfile=args.file_name, threads=args.threads,
                              hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold, blastp_criteria=args.blastp_Evalue_threshold,
                              hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters, blastp_para=args.blastp_parameters)
    else:
        print("Please check the vaild program (genome_table | protein_table | multi_fasta)")

#### main ###
if __name__ == "__main__":
    main()
