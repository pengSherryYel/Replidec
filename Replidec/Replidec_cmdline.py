#!/usr/bin/env python3
# coding: utf-8
# authors: sherry peng, torben sanders, erfan khamespanah
# date: 2026.07.09

import os
import sys
from argparse import RawTextHelpFormatter, ArgumentParser

# Inject parent context directory into sys path to avoid local context resolution issues
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from Replidec.Replidec_multi import bayes_classifier_batch, bayes_classifier_contig, bayes_classifier_genomes
from Replidec.Replidec import checkdb_and_download
from Replidec import __version__

def main():
    """Main execution entry point parsing user flags and routing calls."""
    parser = ArgumentParser(
        description="Replidec: Replication cycle prediction tool for prokaryotic viruses",
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action='version', version=f'Replidec v{__version__}')

    # Fixes an earlier issue by enabling required=False so that database maintenance runs freely without an input file
    parser.add_argument("-p", "--program", default='multi_fasta', required=False,
                        choices=['multi_fasta', 'genome_table', 'protein_table'],
                        metavar="",
                        help=(
                            "{ multi_fasta | genome_table | protein_table }\n\n"
                            "multi_fasta mode:\n"
                            "  input is a fasta file and treats each sequence as one virus\n\n"
                            "genome_table mode:\n"
                            "  input is a tab-separated file [1st col: sample name, 2nd col: genome fasta path]\n\n"
                            "protein_table mode:\n"
                            "  input is a tab-separated file [1st col: sample name, 2nd col: protein faa path]\n"
                        ))

    parser.add_argument("-i", "--input_file", metavar="", default=None,
                        help="The input file, which can be a sequence file or an index table\n")

    parser.add_argument("-w", "--work_dir", default="replidec_results", metavar="", dest="working_directory",
                        help="Directory to store intermediate and final results (default = ./replidec_results)")

    parser.add_argument("-n", "--file_name", default="prediction_summary.tsv", metavar="", dest="file_name",
                        help="Name of final summary file (default = prediction_summary.tsv)")

    parser.add_argument("-t", "--threads", default=4, type=int, metavar="",
                        help="Number of compute threads assigned to execution (default = 4)")

    parser.add_argument("-e", "--hmmer_Eval", default=1e-5, type=float, metavar="", dest="hmmer_Evalue_threshold",
                        help="E-value threshold to filter hmmer results (default = 1e-5)")

    parser.add_argument("-E", "--hmmer_parameters", default="--noali", metavar="", dest="hmmer_parameters",
                        help="Parameters used for hmmer execution (default = --noali)")

    parser.add_argument("-m", "--mmseq_Eval", default=1e-5, type=float, metavar="", dest="mmseqs_Evalue_threshold",
                        help="E-value threshold to filter mmseqs2 results (default = 1e-5)")

    parser.add_argument("-M", "--mmseq_parameters", dest="mmseqs_parameters", metavar="",
                        default="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1",
                        help="Parameters used for mmseqs2 alignment processes")

    parser.add_argument("-b", "--blastp_Eval", default=1e-5, type=float, metavar="", dest="blastp_Evalue_threshold",
                        help="E-value threshold to filter blastp results (default = 1e-5)")

    parser.add_argument("-B", "--blastp_parameter", default="", metavar="", dest="blastp_parameters",
                        help="Parameters used for blastp alignment processes")

    parser.add_argument("-d", "--db_redownload", action='store_true', default=False, dest="db_redownload",
                        help="(Re-)download reference database")

    args = parser.parse_args()
    fileDir = os.path.dirname(os.path.abspath(__file__))

    # Process immediate database installation if requested
    if args.db_redownload:
        print("[INFO] Initiating direct database deployment and verification routines...")
        checkdb_and_download(fileDir, redownload=True)
        if not args.input_file:
            print("[SUCCESS] Database routine finalized successfully. Exiting cleanly.")
            sys.exit(0)

    # Validate that an input file is provided for analysis runs
    if not args.input_file:
        parser.error(
            "the following arguments are required: -i/--input_file (unless only maintaining database environments via -d)")

    print(f"[INFO] Initializing processing pipelines utilizing execution framework: {args.program}")

    # Route execution based on chosen program mode
    if args.program == "genome_table":
        bayes_classifier_genomes(args.input_file, args.working_directory,
                                 summaryfile=args.file_name, threads=args.threads,
                                 hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold,
                                 blastp_criteria=args.blastp_Evalue_threshold,
                                 hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters,
                                 blastp_para=args.blastp_parameters)

    elif args.program == "multi_fasta":
        bayes_classifier_contig(args.input_file, args.working_directory,
                                summaryfile=args.file_name, threads=args.threads,
                                hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold,
                                blastp_criteria=args.blastp_Evalue_threshold,
                                hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters,
                                blastp_para=args.blastp_parameters)

    elif args.program == "protein_table":
        bayes_classifier_batch(args.input_file, args.working_directory,
                               summaryfile=args.file_name, threads=args.threads,
                               hmm_criteria=args.hmmer_Evalue_threshold, mmseqs_criteria=args.mmseqs_Evalue_threshold,
                               blastp_criteria=args.blastp_Evalue_threshold,
                               hmmer_para=args.hmmer_parameters, mmseqs_para=args.mmseqs_parameters,
                               blastp_para=args.blastp_parameters)


if __name__ == "__main__":
    main()
