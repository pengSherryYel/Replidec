#!/usr/bin/env python3
# coding: utf-8
# authors: sherry peng, torben sanders, erfan khamespanah
# date: 2026.07.09

import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

# Inherit the core alignment execution stations, timestamping engines, and global database manifests
from Replidec.Replidec import (
    DATABASE_MANIFEST,
    _get_timestamp,
    checkdb_and_download,
    runProdigal,
    runHmmsearch,
    runMmseqsEasysearch,
    runBlastpsearch,
    load_scoreD,
    load_hmmsearch_opt_batched,
    load_m8_fmt_opt,
    bayes_classifier_single
)
from Replidec.utility import mkdirs


def bayes_classifier_contig(inputfile, wd, summaryfile="BC_predict.summary", threads=4,
                            hmm_criteria=1e-5, mmseqs_criteria=1e-5, blastp_criteria=1e-5,
                            blastp_para="", hmmer_para="", mmseqs_para=""):
    '''
    Aim: High-speed single-pass batch prediction for multi-FASTA contig files.
    '''
    print(f"{_get_timestamp()} [INFO] Initializing reference database verification...")
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)

    if not inputfile:
        print(f"{_get_timestamp()} [CRITICAL] No input file provided.")
        return
    if not os.path.exists(inputfile):
        raise FileNotFoundError(f"Input file '{inputfile}' does not exist.")

    mkdirs(wd)

    # --- Thread Injection Logic ---
    # We dynamically intercept or append the master --threads configuration flag
    # into the underlying string parameters to guarantee underlying binary compliance.
    if '--threads' in mmseqs_para:
        mmseqs_para = re.sub(r'--threads\s+\d+', f'--threads {threads}', mmseqs_para)
    else:
        mmseqs_para += f' --threads {threads}'

    if '--cpu' in hmmer_para:
        hmmer_para = re.sub(r'--cpu\s+\d+', f'--cpu {threads}', hmmer_para)
    else:
        hmmer_para += f' --cpu {threads}'

    if '-num_threads' in blastp_para:
        blastp_para = re.sub(r'-num_threads\s+\d+', f'-num_threads {threads}', blastp_para)
    else:
        blastp_para += f' -num_threads {threads}'

    # --- Phase 1: Metagenomic Gene Prediction ---
    print(f"{_get_timestamp()} [INFO] Sanitizing multi-FASTA boundaries and running Prodigal...")
    sanitized_fasta = os.path.join(wd, "all_contigs_sanitized.fna")
    contig_lengths = {}

    # Standardize header IDs: replace pipeline breaking pipe characters '|' with clean underscores
    with open(sanitized_fasta, "w") as out_f:
        for record in SeqIO.parse(inputfile, "fasta"):
            clean_id = record.id.replace("|", "_")
            record.id = clean_id
            record.description = ""
            contig_lengths[clean_id] = len(record.seq)
            SeqIO.write(record, out_f, "fasta")

    # Run Prodigal in Meta mode to fast-track parallel ORF gene calls across all contigs
    prodigal_wd = os.path.join(wd, "BC_prodigal")
    master_faa = runProdigal(sanitized_fasta, "all_contigs", prodigal_wd, program="meta", otherPara="-g 11")

    if not os.path.exists(master_faa) or os.path.getsize(master_faa) == 0:
        print(f"{_get_timestamp()} [WARNING] No proteins predicted. Exiting workflow cleanly.")
        return

    # Create mapping matrix correlating every predicted ORF back to its parental source contig
    protein_to_contig = {}
    contig_to_proteins = {cid: [] for cid in contig_lengths.keys()}
    for record in SeqIO.parse(master_faa, "fasta"):
        last_underscore = record.id.rfind("_")
        if last_underscore != -1:
            contig_id = record.id[:last_underscore]
            if contig_id in contig_to_proteins:
                protein_to_contig[record.id] = contig_id
                contig_to_proteins[contig_id].append(record.id)

    # --- Phase 2: Vectorized Alignment against Structural Profiles ---
    print(f"{_get_timestamp()} [INFO] Launching vectorized HMMER and MMseqs2 structural profiles...")
    pfam_wd = os.path.join(wd, "BC_pfam")
    integrase_hmm = os.path.join(fileDir, DATABASE_MANIFEST["integrase_hmm"])
    excisionase_hmm = os.path.join(fileDir, DATABASE_MANIFEST["excisionase_hmm"])

    inte_opt = runHmmsearch(master_faa, "all_contigs.BC_integrase", pfam_wd, integrase_hmm, hmmer_para)
    excision_opt = runHmmsearch(master_faa, "all_contigs.BC_excisionase", pfam_wd, excisionase_hmm, hmmer_para)

    # Process alignment results straight into high-speed memory lookups (O(1) lookups)
    inte_anno = load_hmmsearch_opt_batched(inte_opt, criteria=hmm_criteria)
    excision_anno = load_hmmsearch_opt_batched(excision_opt, criteria=hmm_criteria)

    bc_mmseqsDB = os.path.join(fileDir, DATABASE_MANIFEST["mmseqs_index"])
    mmseqs_wd = os.path.join(wd, "BC_mmseqs")
    mmseq_opt = runMmseqsEasysearch(master_faa, "all_contigs.BC_mmseqs", mmseqs_wd, bc_mmseqsDB, otherPara=mmseqs_para)
    mmseq_hits = load_m8_fmt_opt(mmseq_opt, criteria=mmseqs_criteria)

    score_file = os.path.join(fileDir, DATABASE_MANIFEST["scoring_matrix"])
    score_matrix = load_scoreD(score_file)
    p_prior_temperate, p_prior_lytic = score_matrix.get("Prior_probability", [0.0, 0.0])

    # --- Phase 3: Vectorized Chronic Inovirus Marker Screening ---
    print(f"{_get_timestamp()} [INFO] Screening consolidated protein catalog for Chronic Inovirus elements...")
    inno_wd = os.path.join(wd, "BC_Inno")
    inno_hmmDB = os.path.join(fileDir, DATABASE_MANIFEST["inovirus_hmm"])
    inno_blastPre = os.path.join(fileDir, DATABASE_MANIFEST["inovirus_blast"])

    inno_blast_opt = runBlastpsearch(master_faa, "all_contigs.Inno", inno_wd, inno_blastPre,
                                     evalue=str(blastp_criteria), otherPara=blastp_para)
    inno_blast_hits = load_m8_fmt_opt(inno_blast_opt, float(blastp_criteria))

    inno_hmm_opt = runHmmsearch(master_faa, "all_contigs.Inno", inno_wd, inno_hmmDB, hmmer_para)
    inno_hmm_hits = load_hmmsearch_opt_batched(inno_hmm_opt, float(blastp_criteria))

    # --- Phase 4: Scoring Matrix Assembly & Final Report Output ---
    print(f"{_get_timestamp()} [INFO] Assembling scoring matrices and writing global execution summary...")
    summary_path = os.path.join(wd, summaryfile)

    total_contigs = len(contig_lengths)
    processed_count = 0

    with open(summary_path, "w") as opt:
        header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\n"
        opt.write(header)

        for contig_id in contig_lengths.keys():
            proteins = contig_to_proteins[contig_id]

            # Assess core structural markers
            inte_count = sum(1 for p in proteins if p in inte_anno)
            excision_count = sum(1 for p in proteins if p in excision_anno)
            pfam_label = "Temperate" if (inte_count > 0 or excision_count > 0) else "Virulent"

            p_temperate = [p_prior_temperate]
            p_lytic = [p_prior_lytic]
            match_gene_number = 0

            # Accumulate Naive Bayes scores for predicted proteins
            for p in proteins:
                if p in mmseq_hits:
                    match_gene_number += 1
                    ref, _ = mmseq_hits[p]
                    if ref in score_matrix:
                        pt, pl = score_matrix[ref]
                        if float(pt) != 0 and float(pl) != 0:
                            p_temperate.append(pt)
                            p_lytic.append(pl)

            # Compute log-likelihood labels
            if len(p_temperate) == 1:
                # If no proteins match the db
                p_total_temperate, p_total_lytic = 0.0, 0.0
                bc_label = "Unclassified"
            else:
                p_total_temperate = sum(float(t) for t in p_temperate)
                p_total_lytic = sum(float(l) for l in p_lytic)
                if p_total_temperate > p_total_lytic:
                    bc_label = "Temperate"
                elif p_total_temperate < p_total_lytic:
                    bc_label = "Virulent"
                else:
                    bc_label = "Unclassified"

            # Resolve overlapping prediction properties
            if pfam_label == "Virulent" and bc_label == "Unclassified":
                final_label = "Unclassified"
            elif pfam_label == "Temperate" or bc_label == "Temperate":
                final_label = "Temperate"
            else:
                final_label = "Virulent"

            # Overlay Chronic lifestyle flag if Inovirus indicators are caught
            has_inno = any(p in inno_blast_hits or p in inno_hmm_hits for p in proteins)
            if has_inno:
                final_label = "Chronic"

            row = [contig_id, inte_count, excision_count, pfam_label, p_total_temperate, p_total_lytic, bc_label,
                   final_label, match_gene_number]
            opt.write("\t".join(str(x) for x in row) + "\n")

            processed_count += 1
            if processed_count % 5000 == 0 or processed_count == total_contigs:
                print(f"{_get_timestamp()} [PROGRESS] Evaluated {processed_count:,} / {total_contigs:,} contigs.")

    print(f"{_get_timestamp()} [SUCCESS] Job finalized cleanly. Summary written to: {summary_path}")


def bayes_classifier_batch(inputfile, wd, summaryfile="BC_predict.summary", threads=10,
                           hmm_criteria=1e-5, mmseqs_criteria=1e-5, blastp_criteria=1e-5,
                           blastp_para="", hmmer_para="", mmseqs_para=""):
    '''
    Aim: Multithreaded engine execution for index data tables mapping pre-computed protein FAA listings.
    '''
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)
    mkdirs(wd)

    with open(os.path.join(wd, summaryfile), "w") as opt:
        header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
        opt.write(header)

        executor = ThreadPoolExecutor(max_workers=threads)
        all_task = []
        kwargsD = {"hmm_criteria": hmm_criteria, "mmseqs_criteria": mmseqs_criteria, "blastp_criteria": blastp_criteria,
                   "blastp_para": blastp_para, "hmmer_para": hmmer_para, "mmseqs_para": mmseqs_para}

        with open(inputfile) as f:
            for line in f:
                sample_name, path = line.strip("\n").split("\t")
                # REMOVED: time.sleep(10) bottleneck which artificially slowed down multi-sample benchmarks
                all_task.append(executor.submit(bayes_classifier_single, path, sample_name, wd, **kwargsD))

            for future in as_completed(all_task):
                res = future.result()
                res.append(path)
                opt.write("\t".join([str(i) for i in res]) + "\n")


def bayes_classifier_genomes(inputfile, wd, summaryfile="BC_predict.summary", threads=10,
                             hmm_criteria=1e-5, mmseqs_criteria=1e-5, blastp_criteria=1e-5,
                             blastp_para="", hmmer_para="", mmseqs_para=""):
    '''
    Aim: Multithreaded engine execution for structured tables mapping whole genome assembly configurations.
    '''
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)
    mkdirs(wd)

    with open(os.path.join(wd, summaryfile), "w") as opt:
        header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
        opt.write(header)

        executor = ThreadPoolExecutor(max_workers=threads)
        all_task = []
        kwargsD = {"hmm_criteria": hmm_criteria, "mmseqs_criteria": mmseqs_criteria, "blastp_criteria": blastp_criteria,
                   "blastp_para": blastp_para, "hmmer_para": hmmer_para, "mmseqs_para": mmseqs_para}
        faaDict = {}

        with open(inputfile) as f:
            for line in f:
                sample_name, path = line.strip("\n").split("\t")
                # Perform gene calls independently per genome workspace
                faaFile = runProdigal(path, sample_name, f"{wd}/BC_prodigal", program="meta", otherPara="-g 11")
                faaDict[sample_name] = faaFile

            for sample_name, faaFile in faaDict.items():
                if os.path.getsize(faaFile) != 0:
                    # REMOVED: time.sleep(10) bottleneck which introduced dead stalls during task parsing
                    all_task.append(executor.submit(bayes_classifier_single, faaFile, sample_name, wd, **kwargsD))

            for future in as_completed(all_task):
                res = future.result()
                faaFile = faaDict[res[0]]
                res.append(faaFile)
                opt.write("\t".join([str(i) for i in res]) + "\n")


if __name__ == "__main__":
    pass