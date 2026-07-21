#!/usr/bin/env python3
# coding: utf-8
# authors: sherry peng, torben sanders, erfan khamespanah
# mail: xue.peng@helmholtz-muenchen.de
# date: 2026.07.09

import re
import sys
from collections import defaultdict
import os
from subprocess import Popen, PIPE
import time
from Replidec.utility import mkdirs, checkEnv

# Master remote references pointing to static versioned database checkpoints
DATABASE_MANIFEST = {
    "version": "0.3.2",
    "url": "https://zenodo.org/records/15781219/files/db_v0.3.2.tar.gz",
    "integrase_hmm": "db/integrase_pfv34.hmm",
    "excisionase_hmm": "db/excisionase_pfv34.hmm",
    "mmseqs_index": "db/bayes_mmseqs_index/training_prot_04_2025",
    "scoring_matrix": "db/prokaryote_only_training_cluster_04_2025.stat.scoreOpt.tsv",
    "inovirus_hmm": "db/Final_marker_morph.hmm",
    "inovirus_blast": "db/Marker_ALV1"
}

def _get_timestamp():
    """Generates standard tracking prefixes for log reporting."""
    return time.strftime("[%Y-%m-%d %H:%M:%S]")


def runProdigal(inputseq, prefix, wd, program="meta", otherPara="-g 11"):
    """
    Runs Prodigal for gene calling.
    Redirects stdout/stderr to an isolated log within the user's working directory.
    """
    checkEnv("prodigal")
    mkdirs(wd)
    log_file = os.path.join(wd, "external_tools.log")

    cmd = "prodigal -i {0} {4} -a {2}/{1}.prodigal.gene.faa -d {2}/{1}.prodigal.gene.ffn -p {3} -f gff -o {2}/{1}.temp >> {5} 2>&1".format(
        inputseq, prefix, wd, program, otherPara, log_file
    )

    print(f"{_get_timestamp()} [INFO] Running Prodigal: {prefix}")
    obj = Popen(cmd, shell=True)
    obj.wait()
    return os.path.join(wd, f"{prefix}.prodigal.gene.faa")


def runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara="--noali --cpu 1"):
    """Executes target profile scans against reference hidden Markov databases."""
    checkEnv("hmmsearch")
    mkdirs(wd)
    log_file = os.path.join(wd, "external_tools.log")
    output_tbl = os.path.join(wd, f"{prefix}.hmmsearch.tblout")
    output_out = os.path.join(wd, f"{prefix}.hmmsearch.out")

    cmd = "hmmsearch {4} -o {5} --tblout {2} {3} {0} >> {1} 2>&1".format(
        inputfile, log_file, output_tbl, hmmModel, otherPara, output_out
    )

    print(f"{_get_timestamp()} [INFO] Running HMMer search: {prefix}")
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output_tbl


def runMmseqsEasysearch(inputfile, prefix, wd, hmmDB,
                        otherPara="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 "
                                  "--min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 30"):
    """Executes fast protein homology searches using MMseqs2."""
    checkEnv("mmseqs")
    mkdirs(wd)
    log_file = os.path.join(wd, "external_tools.log")

    tmpdir = os.path.join(wd, f"{prefix}_tmp")
    mkdirs(tmpdir)
    output = os.path.join(wd, f"{prefix}.mmseqs.m8")

    cmd = "mmseqs easy-search {0} {1} {2} {3} {4} >> {5} 2>&1 && rm -rf {3}".format(
        inputfile, hmmDB, output, tmpdir, otherPara, log_file
    )

    print(f"{_get_timestamp()} [INFO] Running MMseqs2 easy-search: {prefix}")
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def runBlastpsearch(inputfile, prefix, wd, db_prefix, evalue="1e-3", otherPara="-num_threads 3"):
    """Standard BLASTp sequence alignment."""
    checkEnv("blastp")
    mkdirs(wd)
    log_file = os.path.join(wd, "external_tools.log")
    output = os.path.join(wd, f"{prefix}.blastp.m8")

    cmd = "blastp -query {0} -db {1} -out {2} -outfmt 6 -evalue {3} {4} >> {5} 2>&1".format(
        inputfile, db_prefix, output, evalue, otherPara, log_file
    )

    print(f"{_get_timestamp()} [INFO] Running BLASTp search: {prefix}")
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def load_scoreD(score_file):
    """
    Parses the pre-computed cluster conditional log-probability matrix.
    Returns: Dict mapping a reference element to its [temperate_score, lytic_score].
    """
    d = {}
    with open(score_file) as f:
        for line in f:
            ref, temperate, virulent, members = line.strip("\n").split("\t")
            for member in members.split(";"):
                d[member] = [temperate, virulent]
    return d


def load_hmmsearch_opt(hmmsearch_opt, criteria=1e-5):
    """
    Parses sequential single-sample HMMER --tblout files.
    """
    annoD = defaultdict(dict)
    if not os.path.exists(hmmsearch_opt):
        return annoD
    with open(hmmsearch_opt) as f:
        for line in f:
            if not line.startswith("#"):
                t = re.split(r"\s+", line.strip("\n"))
                if len(t) < 5:
                    continue
                target_name, _, query_name, _, Evalue, *rest = t
                bst_Evalue = rest[2] if len(rest) > 2 else Evalue
                accession = t[3].split(".")[0]

                if float(Evalue) <= float(criteria) and float(bst_Evalue) <= float(criteria):
                    annoD[target_name][accession] = query_name
        return annoD


def load_hmmsearch_opt_batched(hmmsearch_opt, criteria=1e-5):
    """
    High-speed parser for batch/consolidated HMMER tabular outputs.
    Extracts column 1 (target protein sequence identifiers) into an optimized set.
    """
    hit_proteins = set()
    if not os.path.exists(hmmsearch_opt):
        return hit_proteins

    with open(hmmsearch_opt) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = re.split(r"\s+", line.strip("\n"))
            if len(t) < 5:
                continue

            # Column mapping: t[0] yields target sequence ID, t[2] holds query matrix name
            target_name, _, query_name, _, Evalue = t[:5]

            if float(Evalue) <= float(criteria):
                hit_proteins.add(target_name)

    return hit_proteins


def load_m8_fmt_opt(m8_input, criteria=1e-5):
    """
    Parses blastp / mmseqs2 blast-style m8 (-outfmt 6) layout structures.
    Filters by E-value and maps query sequences to their top hit reference IDs.
    """
    d = {}
    if os.path.exists(m8_input):
        with open(m8_input) as f:
            for line in f:
                parts = line.strip("\n").split("\t")
                if len(parts) < 11:
                    continue
                query, ref, _, _, _, _, _, _, _, _, evalue = parts[:11]
                if float(evalue) <= criteria:
                    if query not in d or float(evalue) < float(d[query][-1]):
                        d[query] = [ref, evalue]
    return d


def calculate_score(mmseqOpt, scoreD, criteria=1e-5):
    """
    Computes Naive Bayes joint conditional probabilities for lifestyles.
    Combines log-likelihood values from homological alignments.
    """
    p_prior_temperate, p_prior_lytic = scoreD.get("Prior_probability", [0.0, 0.0])

    p_temperate = [p_prior_temperate]
    p_lytic = [p_prior_lytic]

    d = load_m8_fmt_opt(mmseqOpt, criteria)
    match_gene_number = len(d)

    for query, values in d.items():
        ref, evalue = values
        if ref in scoreD:
            p_temperate_gc, p_lytic_gc = scoreD[ref]
            if float(p_temperate_gc) != 0 and float(p_lytic_gc) != 0:
                p_temperate.append(p_temperate_gc)
                p_lytic.append(p_lytic_gc)

    label = "NA"
    p_total_temperate, p_total_lytic = 0.0, 0.0
    if len(p_temperate) == 1:
        label = "Unclassified"
    else:
        for t, l in zip(p_temperate, p_lytic):
            p_total_temperate += float(t)
            p_total_lytic += float(l)

        if p_total_temperate > p_total_lytic:
            label = "Temperate"
        elif p_total_temperate < p_total_lytic:
            label = "Virulent"
        else:
            label = "Unclassified"

    return p_total_temperate, p_total_lytic, label, match_gene_number


def inoviruses_PI_like_gene_search(inputfile, prefix, wd, hmmdb, blastdb_prefix, hmm_evalue="1e-3",
                                   blastp_evalue="1e-3", blastp_para="-num_threads 3", hmmer_para="--noali --cpu 1"):
    """Detects chronic Inovirus signature sequences using single-sample profiles."""
    inovBlastpOpt = runBlastpsearch(inputfile, prefix, wd, blastdb_prefix, evalue=blastp_evalue, otherPara=blastp_para)
    innoBlastD = load_m8_fmt_opt(inovBlastpOpt, float(blastp_evalue))

    innoHmmerOpt = runHmmsearch(inputfile, prefix, wd, hmmdb, otherPara=hmmer_para)
    innoHmmerD = load_hmmsearch_opt(innoHmmerOpt, float(hmm_evalue))

    return bool(innoBlastD or innoHmmerD)


def check_db_md5(db_dir):
    """Validates the local database files against official MD5 checksum targets."""
    cmd = "cd %s/db && md5sum --check md5sum.list" % (db_dir)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    obj.wait()
    out, err = obj.communicate()
    if b"FAILED" in out or b"FAILED" in err:
        print(f"{_get_timestamp()} [WARNING] Database validation failure detected! Re-download recommended.")
        sys.exit(1)
    else:
        print(f"{_get_timestamp()} [INFO] Database md5 integrity checks: PASSED")


def checkdb_and_download(scriptPos, redownload=False):
    """Ensures reference database components are present and structurally valid."""
    if redownload:
        if os.path.exists(os.path.join(scriptPos, "db")):
            cmd = f"rm -rf {os.path.join(scriptPos, 'discarded_db')} && mv -f {os.path.join(scriptPos, 'db')} {os.path.join(scriptPos, 'discarded_db')}"
            obj = Popen(cmd, shell=True)
            obj.wait()

    if not os.path.exists(os.path.join(scriptPos, "db")):
        print(f"{_get_timestamp()} [INFO] Reference database missing. Downloading dependencies...")
        url = DATABASE_MANIFEST["url"]
        file = os.path.split(url)[-1]
        cmd = "wget {0} -P {1} && cd {1} && tar -zxvf {2} && rm -rf {2}".format(url, scriptPos, file)

        try:
            obj = Popen(cmd, shell=True)
            obj.wait()
            check_db_md5(scriptPos)
        except Exception as e:
            print(f"{_get_timestamp()} [CRITICAL] Database setup failed: {str(e)}. Rerun with -d.")
            sys.exit(1)
    else:
        print(f"{_get_timestamp()} [INFO] Database verified.")


def bayes_classifier_single(inputfile, prefix, wd, hmm_criteria=1e-5, mmseqs_criteria=1e-5, blastp_criteria=1e-3,
                            blastp_para="-num_threads 3", hmmer_para="--noali --cpu 3",
                            mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):
    """
    Processes single-sample lifestyle predictions.
    Retained to support alternative processing workflows.
    """
    mkdirs(wd)
    fileDir = os.path.dirname(os.path.abspath(__file__))

    # --- Step 1: Pre-Screen for Chronic Markers to Short-Circuit Early ---
    inno_hmmDB = os.path.join(fileDir, "db/Final_marker_morph.hmm")
    inno_blastPre = os.path.join(fileDir, "db/Marker_ALV1")
    inno_dect_wd = os.path.join(wd, "BC_Inno")

    inno_res = inoviruses_PI_like_gene_search(
        inputfile, f"{prefix}.Inno", inno_dect_wd, inno_hmmDB, inno_blastPre,
        hmm_evalue=blastp_criteria, blastp_evalue=blastp_criteria, blastp_para=blastp_para, hmmer_para=hmmer_para
    )
    
    if inno_res:
        # Short circuit: return NA's for scores and "Skipped" for all structural labels, assigning "Chronic" as final label
        return [prefix, "NA", "NA", "Skipped", "NA", "NA", "Skipped", "Chronic", "NA"]

    # --- Step 2: Proceed With Structural Alignments if Not Chronic ---
    pfam_label, bc_label, final_label = "Unclassified", "Unclassified", "Unclassified"
    integrase_hmm = os.path.join(fileDir, "db/integrase_pfv34.hmm")
    excisionase_hmm = os.path.join(fileDir, "db/excisionase_pfv34.hmm")

    pfam_wd = os.path.join(wd, "BC_pfam")
    inte_opt = runHmmsearch(inputfile, f"{prefix}.BC_integrase", pfam_wd, integrase_hmm, hmmer_para)
    excision_opt = runHmmsearch(inputfile, f"{prefix}.BC_excisionase", pfam_wd, excisionase_hmm, hmmer_para)

    inte_annoD = load_hmmsearch_opt(inte_opt, criteria=hmm_criteria)
    excision_annoD = load_hmmsearch_opt(excision_opt, criteria=hmm_criteria)

    inte_label = len(inte_annoD) if inte_annoD else 0
    excision_label = len(excision_annoD) if excision_annoD else 0

    if inte_label or excision_label:
        pfam_label = "Temperate"
    else:
        pfam_label = "Virulent"
        
    bc_mmseqsDB = os.path.join(fileDir, "db/bayes_mmseqs_index/training_prot_04_2025")
    mmseqs_wd = os.path.join(wd, "BC_mmseqs")
    mmseq_opt = runMmseqsEasysearch(inputfile, f"{prefix}.BC_mmseqs", mmseqs_wd, bc_mmseqsDB, otherPara=mmseqs_para)

    score_file = os.path.join(fileDir, "db/prokaryote_only_training_cluster_04_2025.stat.scoreOpt.tsv")
    member2scoreD = load_scoreD(score_file)

    p_total_temperate, p_total_lytic, bc_label, match_gene_number = calculate_score(mmseq_opt, member2scoreD,
                                                                                    criteria=mmseqs_criteria)

    if pfam_label == "Temperate" or bc_label == "Temperate":
        final_label = "Temperate"
    elif pfam_label == "Virulent" and bc_label == "Unclassified":
        final_label = "Unclassified"
    else:
        final_label = "Virulent"

    return [prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label,
            match_gene_number]


if __name__ == "__main__":
    scriptPos = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(scriptPos)
