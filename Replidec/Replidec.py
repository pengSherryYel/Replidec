#!/usr/bin/env python3
# coding: utf-8
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6


import re, sys
from collections import defaultdict
import os
from subprocess import Popen, PIPE
import math
from Replidec.utility import mkdirs, checkEnv
from Bio import Seq, SeqIO
from threading import Thread


# ------ function ------
def runProdigal(
        inputseq,
        prefix,
        wd,
        program="meta",
        otherPara="-g 11"
):

    """
    Aim: run prodigal to predict gene from contig

    Usage: runProdigal(inputseq,prefix,wd,program="meta")
        inputfile: input seq file
        prefix: output prefix
        wd: work path where put the result
        program: euqal to -p parameter in prodigal(meta(default)|single)
        otherPara: parameter for run the prodigal
            default: "-g 11"

    Return: output file path (*.faa)
    """

    checkEnv("prodigal")
    mkdirs(wd)
    cmd = "prodigal -i {0} {4} -a {2}/{1}.prodigal.gene.faa -d {2}/{1}.prodigal.gene.ffn -p {3} -f gff -o {2}/{1}.temp 2>&1 >>log".format(
        inputseq, prefix, wd, program, otherPara
    )
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    print("prodigal done!")
    return "%s/%s.prodigal.gene.faa" % (wd, prefix)


def runHmmsearch(
        inputfile,
        prefix,
        wd,
        hmmModel,
        otherPara="--noali --cpu 1"
):

    """
    Aim: run hmmer search for a pfam hmm database

    Usage: runHmmsearch(inputfile,prefix,wd,hmmModel,otherPara="--cpu 1")
        inputfile: a protein set from a metabin or a genome
        prefix: sample identifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmModel: a pfam *.hmm file
        otherPara: parameter for run the hmmer search.
            default: "--cpu 1"

    Return: output file path (*.tblout)
    """

    checkEnv("hmmsearch")
    mkdirs(wd)
    cmd = "hmmsearch {4} -o {2}/{1}.hmmsearch.out --tblout {2}/{1}.hmmsearch.tblout {3} {0} 2>&1 >>log".format(
        inputfile, prefix, wd, hmmModel, otherPara
    )
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    print("hmmsearch done!")
    return "%s/%s.hmmsearch.tblout" % (wd, prefix)


def runMmseqsEasysearch(
    inputfile,
    prefix,
    wd,
    hmmDB,
    otherPara="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 30",
):
    """
    Aim: run mmseq easy search for a fasta database

    Usage: runMmseqsEasysearch(inputfile,prefix,wd,hmmDB,otherPara)
        inputfile: a protein set from a metabin or a genome
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmDB: a fasta database built from mmseqs creatdb command
        otherPara: parameter for run the mmseqs easy search.
            default: "-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40
                      --cov-mode 0 --greedy-best-hits 1 --threads 30"

    Return: output file path (*.m8)
    """

    checkEnv("mmseqs")
    mkdirs(wd)

    tmpdir = "%s_tmp" % prefix
    mkdirs(tmpdir)
    output = os.path.join(wd, "%s.mmseqs.m8" % prefix)

    cmd = "mmseqs easy-search {0} {1} {2} {3} {4} 2>&1 >>log && rm -rf {3} && echo 'mmesqs search done!'".format(
        inputfile, hmmDB, output, tmpdir, otherPara
    )
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def runBlastpsearch(
    inputfile,
    prefix,
    wd,
    db_prefix,
    evalue="1e-3",
    otherPara="-num_threads 3"
):

    """
    Aim: run blastp search for a fasta database(db should build first)

    Usage: runBlastpsearch(inputfile,db,prefix,wd,evalue,otherPara)
        inputfile: a protein set from a metabin or a genome
        db_prefix: database prefix for the blast(db should build before)
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        evalue: evaule used in blastp command (default: 1e-3)
        otherPara: parameter for run the blastp search. (default: -num_threads 3)

    Return: output file path (*.m8)
    """

    checkEnv("blastp")
    mkdirs(wd)

    output = os.path.join(wd, "%s.blastp.m8" % prefix)

    cmd = "blastp -query {0} -db {1} -out {2} -outfmt 6 -evalue {3} {4} 2>&1 >>log".format(
        inputfile, db_prefix, output, evalue, otherPara
    )
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def load_scoreD(score_file):

    """
    Aim: Parse the score file
    Return: dict. d[member] = [temperate_score, lytic_score]
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
    Aim: parse the hmmersearch output. this file contain multiple columns. is the output from -tblout parameter

    Return: dict.  d[refname][annoacc] = description
    """

    print("loading mmsearch output")
    annoD = defaultdict(dict)
    with open(hmmsearch_opt) as f:
        for line in f:
            if not line.startswith("#"):
                t = re.split(r"\s+", line.strip("\n"))
                (
                    target_name,
                    target_accession,
                    query_name,
                    accession,
                    Evalue,
                    score,
                    bias,
                    bst_Evalue,
                    bst_score,
                    bst_bias,
                    exp,
                    reg,
                    clu,
                    ov,
                    env,
                    dom,
                    rep,
                    inc,
                    *description_of_target,
                ) = t
                accession = accession.split(".")[0]

                if float(Evalue) <= float(criteria) and float(bst_Evalue) <= float(
                    criteria
                ):
                    annoD[target_name][accession] = query_name
        return annoD


def load_m8_fmt_opt(m8_input, criteria=1e-5):
    """
    Aim: parse m8 format (output from mmseqs and blastp)

    Return: dict.  d[query] = [ref,evalue]
    """

    d = {}
    if os.path.exists(m8_input):
        with open(m8_input) as f:
            for line in f:
                (
                    query,
                    ref,
                    iden,
                    length,
                    mismatch,
                    gap,
                    qstart,
                    qend,
                    sstart,
                    send,
                    evalue,
                    bit_score,
                ) = line.split("\t")
                if query not in d:
                    if float(evalue) <= criteria:
                        d[query] = [ref, evalue]
                else:
                    if float(evalue) < float(d[query][-1]):
                        d[query] = [ref, evalue]
    return d


def inoviruses_PI_like_gene_search(
    inputfile,
    prefix,
    wd,
    hmmdb,
    blastdb_prefix,
    hmm_evalue="1e-3",
    blastp_evalue="1e-3",
    blastp_para="-num_threads 3",
    hmmer_para="--noali --cpu 1",
):
    """
    Aim: search inoviruses PI like gene using hmmsearch and blastp

    Return: Bool value(TRUE: contain inoviruses marker gene.)
    """

    inovBlastpOpt = runBlastpsearch(
        inputfile,
        prefix,
        wd,
        blastdb_prefix,
        evalue=blastp_evalue,
        otherPara=blastp_para,
    )
    innoBlastD = load_m8_fmt_opt(inovBlastpOpt, blastp_evalue)

    ## run mmseqs
    innoHmmerOpt = runHmmsearch(inputfile, prefix, wd, hmmdb, otherPara=hmmer_para)
    innoHmmerD = load_hmmsearch_opt(innoHmmerOpt, hmm_evalue)

    if innoBlastD or innoHmmerD:
        return True
    else:
        return False


def calcaulate_score(mmseqOpt, scoreD, criteria=1e-5):

    """
    Aim: calculate P(temperate|GC1,GC2,...,GCN) and P(virulent||GC1,GC2,...,GCN) for a given mmseqOpt

    Usage: calculate_score(mmseqOpt,scoreD,criteria=1e-5)
        mmseqOpt: mmseq easy search output with our database from function runMmseqsEasysearch()
        scoreD: the return dict from function load_scoreD()
        criteria: just select e-value greater than criteria

    Return: p_total_temperate,p_total_lytic,label
    """


    p_prior_temperate, p_prior_lytic = scoreD.get("Prior_probability", "NA")

    p_temperate = []
    p_lytic = []
    p_temperate.append(p_prior_temperate)
    p_lytic.append(p_prior_lytic)

    # find the smallest e-value result for the query => return a dict
    d = load_m8_fmt_opt(mmseqOpt, criteria)

    # calculate how many genes are matched
    match_gene_number = len(d)

    # enumerate the dict to add the score of each lifestyle to lists
    for query, values in d.items():
        ref, evalue = values
        p_temperate_gc, p_lytic_gc = scoreD[ref]
        if p_temperate_gc != 0 and p_lytic_gc != 0:
            p_temperate.append(p_temperate_gc)
            p_lytic.append(p_lytic_gc)

    # calculate P(temperate|GC1,GC2,...,GCN) and P(virulent||GC1,GC2,...,GCN)
    label = "NA"
    p_total_temperate, p_total_lytic = 0, 0
    if len(p_temperate) == 1:
        label = "NotEnoughInfo"
    else:
        for t, l in zip(p_temperate, p_lytic):
            p_total_temperate += float(t)
            p_total_lytic += float(l)

    if p_total_temperate > p_total_lytic:
        label = "Temperate"
    elif p_total_temperate < p_total_lytic:
        label = "Virulent"
    else:
        label = "NotEnoughInfo"
    return p_total_temperate, p_total_lytic, label, match_gene_number


def chunk_list(inputlist, chunksize=10):
    d = {}
    n = 0
    for i in range(0, len(inputlist), chunksize):
        t = inputlist[i : i + chunksize]
        d[n] = t
        n += 1
    return d


def check_db_md5(db_dir):
    cmd = "cd %s/db && md5sum --check md5sum.list && cd -" % (db_dir)
    obj = Popen(cmd, shell=True, stdout=PIPE)
    obj.wait()

    print([i.decode("utf-8") for i in obj.stdout.readlines()])
    l = [i.decode("utf-8") for i in obj.stdout.readlines()]
    if "FAILED" in ";".join([i.strip("\n") for i in l]):
        print("WARNING: Please recheck the databse, database file have wrong information.")
        sys.exit()
    else:
        print("db check done!")


def checkdb_and_download(scriptPos, redownload=False):
    if redownload:
        mkdirs("discarded_db")
        if os.path.exists("%s/db" % scriptPos):
            cmd = "mv -f %s/db discarded_db " % scriptPos
            obj = Popen(cmd, shell=True)
            obj.wait()

    if not os.path.exists(os.path.join(scriptPos, "db")):
        print("db not exist! download database ...")
        url = "https://zenodo.org/records/15781219/files/db_v0.3.2.tar.gz"
        file = os.path.split(url)[-1]
        cmd = "wget {0} -P {1} && cd {1} && tar -zxvf {2} && rm -rf {2}".format(
            url, scriptPos, file
        )
        
        ## add error control for db download
        try:
            obj = Popen(cmd, shell=True, stdout=PIPE)
            obj.wait()

            check_db_md5(scriptPos)
        except:
            print("db download error!! Please RERUN with -d")
            sys.exit()
        else:
            print("db download done!")

    else:
        print("db exist!")


def bayes_classifier_single(
    inputfile,
    prefix,
    wd,
    hmm_criteria=1e-5,
    mmseqs_criteria=1e-5,
    blastp_criteria=1e-3,
    blastp_para="-num_threads 3",
    hmmer_para="--noali --cpu 3",
    mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 3",
):
    """
    Aim: single predict lifestyle

    Process: have 3 phase
        phase 1: align to the integrase and excisionase(16 pfam id in total) -- hmmer
        phase 2: align to our protein database -- mmseqs easy-search
        phase 3: dected innovirues -- hmmer and blastp
        if temperate appear in either one phase --> final will be temperate

    Usage: bayes_classifier_single(inputfile,prefix,wd,hmm_criteria=1e-5,mmseqs_criteria=1e-5)
        inputfile: protein set from a bin|a genome.
        prefix: sample identifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmm_criteria: criteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_criteria: criteria to filter mmseqs evalue greater than x (default: 1e-5)
        blastp_criteria: criteria to filter blastp evalue greater than x (default: 1e-5)

    Return:[prefix,inte_label,excision_label,pfam_label,p_total_temperate,p_total_lytic,bc_label,final_label]
    """
    print("## %s start!" % prefix)
    mkdirs(wd)
    pfam_label, bc_label, final_label = "Virulent", "Virulent", "Virulent"
    fileDir = os.path.dirname(os.path.abspath(__file__))

    # phase 1 hmmsearch for pfam
    integrase_hmm = os.path.join(fileDir, "db/integrase_pfv34.hmm")
    excisionase_hmm = os.path.join(fileDir, "db/excisionase_pfv34.hmm")

    # run hmmsearch for integrease and excisionase
    pfam_wd = os.path.join(wd, "BC_pfam")
    mkdirs(pfam_wd)
    inte_opt = runHmmsearch(
        inputfile, "%s.BC_integrase" % prefix, pfam_wd, integrase_hmm, hmmer_para
    )
    excision_opt = runHmmsearch(
        inputfile, "%s.BC_excisionase" % prefix, pfam_wd, excisionase_hmm, hmmer_para
    )

    # load pfam result
    inte_label, excision_label = [0, 0]
    inte_annoD = load_hmmsearch_opt(inte_opt, criteria=hmm_criteria)
    excision_annoD = load_hmmsearch_opt(excision_opt, criteria=hmm_criteria)

    if inte_annoD:
        inte_label = len(inte_annoD)
    if excision_annoD:
        excision_label = len(excision_annoD)

    if inte_label or excision_label:
        pfam_label = "Temperate"

    # phase 2 bayes classifier
    # run mmseqs search
    print("Using prokaryote only protein v04_2025 as DB")
    bc_mmseqsDB = os.path.join(
            fileDir, "db/bayes_mmseqs_index/training_prot_04_2025"
    )

    mmseqs_wd = os.path.join(wd, "BC_mmseqs")
    mmseqs_prefix = "%s.BC_mmseqs" % prefix
    mmseq_opt = runMmseqsEasysearch(
        inputfile, mmseqs_prefix, mmseqs_wd, bc_mmseqsDB, otherPara=mmseqs_para
    )


    # load score file
    score_file = os.path.join(
            fileDir, "db/prokaryote_only_training_cluster_04_2025.stat.scoreOpt.tsv"
    )

    member2scoreD = load_scoreD(score_file)

    # read mmseqs easy search file
    p_total_temperate, p_total_lytic, bc_label, match_gene_number = calcaulate_score(
        mmseq_opt, member2scoreD, criteria=mmseqs_criteria
    )
    print(
        prefix,
        inte_label,
        excision_label,
        pfam_label,
        p_total_temperate,
        p_total_lytic,
        bc_label,
    )

    # phase 2B combine result pfam and bayes classifier
    if pfam_label == "Temperate" or bc_label == "Temperate":
        final_label = "Temperate"

    # phase 3: detect Innovirues
    # way1: based on marker gene
    inno_hmmDB = os.path.join(fileDir, "db/Final_marker_morph.hmm")
    inno_blastPre = os.path.join(fileDir, "db/Marker_ALV1")
    inno_dect_wd = os.path.join(wd, "BC_Inno")
    inno_prefix = "%s.Inno" % prefix

    ### in this 1e-3 is used to identify the PI-like
    inno_res = inoviruses_PI_like_gene_search(
        inputfile,
        inno_prefix,
        inno_dect_wd,
        inno_hmmDB,
        inno_blastPre,
        hmm_evalue=blastp_criteria,
        blastp_evalue=blastp_criteria,
        blastp_para=blastp_para,
        hmmer_para=hmmer_para,
    )
    if inno_res == True:
        final_label = "Chronic"
    print("## %s end!" % prefix)
    return [
        prefix,
        inte_label,
        excision_label,
        pfam_label,
        p_total_temperate,
        p_total_lytic,
        bc_label,
        final_label,
        match_gene_number,
    ]


if __name__ == "__main__":
    scriptPos = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(scriptPos)
