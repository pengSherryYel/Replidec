#!/usr/bin/env python3
# coding: utf-8
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6


import re,sys
from collections import defaultdict
import os
from subprocess import Popen,PIPE
import math
from Replidec.utility import mkdirs, checkEnv
from Bio import Seq,SeqIO
from threading import Thread

# ------ function ------
def runProdigal(inputseq, prefix, wd, program="meta", otherPara="-g 11"):
    '''
    Aim: run prodigal to predict gene from contig

    Usage: runProdigal(inputseq,prefix,wd,program="meta")
        inputfile: input seq file
	prefix: output prefix
        wd: work path where put the result
	program: euqal to -p parameter in prodigal(meta(default)|single)
        otherPara: parameter for run the prodigal
            default: "-g 11"

    Return: output file path (*.faa)
    '''
    #cmd = "prodigal -a {2}/{1}.prodigal.gene.faa -d {2}/{1}.prodigal.gene.fna -g 11  -i {0} -o {2}/{1}.prodigal.output -s {2}/{1}.prodigal.gene.score".format(inputseq,prefix,wd)
    checkEnv("prodigal")
    mkdirs(wd)
    cmd = "prodigal -i {0} {4} -a {2}/{1}.prodigal.gene.faa -d {2}/{1}.prodigal.gene.ffn -p {3} -f gff -o {2}/{1}.temp 2>&1 >>log".format(
            inputseq, prefix, wd, program, otherPara)
    print("RUN command: %s\n"%cmd)
    obj = Popen(cmd,shell=True)
    obj.wait()
    print("prodigal done!")
    return "%s/%s.prodigal.gene.faa"%(wd,prefix)


def runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara="--noali --cpu 1"):
    '''
    Aim: run hmmer search for a pfam hmm database

    Usage: runHmmsearch(inputfile,prefix,wd,hmmModel,otherPara="--cpu 1")
        inputfile: a protein set from a metabin or a genome
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmModel: a pfam *.hmm file
        otherPara: parameter for run the hmmer search.
            default: "--cpu 1"

    Return: output file path (*.tblout)
    '''

    checkEnv("hmmsearch")
    mkdirs(wd)
    cmd = "hmmsearch {4} -o {2}/{1}.hmmsearch.out --tblout {2}/{1}.hmmsearch.tblout {3} {0} 2>&1 >>log".format(
        inputfile, prefix, wd, hmmModel, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    print("hmmsearch done!")
    return "%s/%s.hmmsearch.tblout" % (wd, prefix)


def runMmseqsEasysearch(inputfile, prefix, wd, hmmDB, otherPara="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 30"):
    '''
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
    '''
    checkEnv("mmseqs")
    mkdirs(wd)

    tmpdir = "%s_tmp" % prefix
    mkdirs(tmpdir)
    output = os.path.join(wd, "%s.mmseqs.m8" % prefix)

    cmd = "mmseqs easy-search {0} {1} {2} {3} {4} 2>&1 >>log && rm -rf {3} && echo 'mmesqs search done!'".format(
        inputfile, hmmDB, output, tmpdir, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output

def runBlastpsearch(inputfile,prefix,wd,db_prefix,evalue="1e-3",otherPara="-num_threads 3"):
    '''
    Aim: run blastp search for a fasta database(db should build first)

    Usage: runBlastpsearch(inputfile,db,prefix,wd,evalue,otherPara)
        inputfile: a protein set from a metabin or a genome
        db_prefix: database prefix for the blast(db should build before)
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        evalue: evaule used in blastp command (default: 1e-3)
        otherPara: parameter for run the blastp search. (default: -num_threads 3)

    Return: output file path (*.m8)
    '''

    checkEnv("blastp")
    mkdirs(wd)

    output = os.path.join(wd, "%s.blastp.m8" % prefix)

    cmd = "blastp -query {0} -db {1} -out {2} -outfmt 6 -evalue {3} {4} 2>&1 >>log".format(
            inputfile,db_prefix,output,evalue,otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def load_scoreD(score_file):
    '''
    Aim: Parse the score file
    Return: dict. d[member] = [temperate_score, lytic_score]
    '''
    d = {}
    with open(score_file) as f:
        for line in f:
            ref, tmperate, virulent, members = line.strip("\n").split("\t")
            for member in members.split(";"):
                d[member] = [tmperate, virulent]
    return d


# def load_mmseq2_cluster_tsv(mmseq2_tsv):
#    d = defaultdict(list)
#    with open(mmseq2_tsv) as f:
#        for line in f:
#            pre,member = line.strip("\n").split("\t")
#            d[pre].append(member)
#    return d


def load_hmmsearch_opt(hmmsearch_opt, creteria=1e-5):
    '''
    Aim: parse the hmmersearch output. this file contain multiple columns. is the output from -tblout parameter

    Return: dict.  d[refname][annoacc] = description
    '''
    print("loading mmsearch output")
    annoD = defaultdict(dict)
    with open(hmmsearch_opt) as f:
        for line in f:
            if not line.startswith("#"):
                t = re.split(r"\s+", line.strip("\n"))
                target_name, target_accession, query_name, accession, Evalue, score, bias, bst_Evalue, bst_score, bst_bias,\
                    exp, reg, clu, ov, env, dom, rep, inc, *description_of_target = t
                accession = accession.split(".")[0]
                # print(target_name,Evalue,bst_Evalue)
                if float(Evalue) <= float(creteria) and float(bst_Evalue) <= float(creteria):
                    annoD[target_name][accession] = query_name
        return annoD

def load_m8_fmt_opt(m8_input, creteria=1e-5):
    '''
    Aim: parse m8 format (output from mmseqs and blastp)

    Return: dict.  d[query] = [ref,evalue]
    '''
    # find the smallest e-value result for the query => return a dict
    d = {}
    if os.path.exists(m8_input):
        with open(m8_input) as f:
            for line in f:
                query, ref, iden, length, mismatch, gap, qstart, qend, sstart, send, evalue, bit_score = line.split("\t")
                if query not in d:
                    if float(evalue) <= creteria:
                        d[query] = [ref, evalue]
                else:
                    if float(evalue) < float(d[query][-1]):
                        d[query] = [ref, evalue]
    return d

def inoviruses_PI_like_gene_search(inputfile,prefix,wd,
        hmmdb,blastdb_prefix,
        hmm_evalue="1e-3", blastp_evalue="1e-3",
        blastp_para="-num_threads 3",
        hmmer_para="--noali --cpu 1"):
    '''
    Aim: search inoviruses PI like gene using hmmsearch and blastp

    Return: Bool value(TRUE: contain inoviruses marker gene.)
    '''
    ## run blastp
    inovBlastpOpt = runBlastpsearch(inputfile,prefix,wd,blastdb_prefix,evalue=blastp_evalue,otherPara=blastp_para)
    innoBlastD = load_m8_fmt_opt(inovBlastpOpt, blastp_evalue)

    ## run mmseqs
    innoHmmerOpt = runHmmsearch(inputfile, prefix, wd, hmmdb, otherPara=hmmer_para)
    innoHmmerD = load_hmmsearch_opt(innoHmmerOpt, hmm_evalue)

    if innoBlastD or innoHmmerD:
        return True
    else:
        return False


def calcaulate_score(mmseqOpt, scoreD, creteria=1e-5):
    '''
    Aim: calculate P(temperate|GC1,GC2,...,GCN) and P(virulent||GC1,GC2,...,GCN) for a given mmseqOpt

    Usage: calcaulate_score(mmseqOpt,scoreD,creteria=1e-5)
        mmseqOpt: mmseq easy seach output with our database from function runMmseqsEasysearch()
        scoreD: the return dict from function load_scoreD()
        creteria: just select e-value greater than creteria

    Return: p_total_temperate,p_total_lytic,label
    '''
    # load the prior probability
    p_prior_temperate, p_prior_lytic = scoreD.get("Prior_probability", "NA")

    p_temperate = []
    p_lytic = []
    p_temperate.append(p_prior_temperate)
    p_lytic.append(p_prior_lytic)

    # find the smallest e-value result for the query => return a dict
    d = load_m8_fmt_opt(mmseqOpt, creteria)
#    d = {}
#    with open(mmseqOpt) as f:
#        for line in f:
#            query, ref, iden, length, mismatch, gap, qstart, qend, sstart, send, evalue, bit_score = line.split(
#                "\t")
#            if query not in d:
#                if float(evalue) <= creteria:
#                    d[query] = [ref, evalue]
#            else:
#                if float(evalue) < float(d[query][-1]):
#                    d[query] = [ref, evalue]

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
        t = inputlist[i:i+chunksize]
        d[n] = t
        n += 1
    return d


#prophage_relate = ['PF00552', 'PF00589', 'PF00665', 'PF02022', 'PF02899',
#                   'PF02920', 'PF09003', 'PF12482', 'PF12834', 'PF12835', 'PF13009', 'PF13102',
#                   'PF13333', 'PF13495', 'PF13683', 'PF13976', 'PF14659', 'PF14882', 'PF16795',
#                   'PF17921', 'PF18103', 'PF18644', 'PF18697', 'PF06806', 'PF07825', 'PF09035']

def check_db_md5(db_dir):
    cmd = "cd %s/db && md5sum --check md5sum.list && cd -"%(db_dir)
    obj = Popen(cmd,shell=True, stdout=PIPE)
    obj.wait()
    print([i.decode("utf-8") for i in obj.stdout.readlines()])
    if "FAILED" in ";".join([i.decode("utf-8") for i in obj.stdout.readlines()]):
        print("Please recheck the databse, database file have wrong information.")
        sys.exit()
    else:
        print("db download done!")

def checkdb_and_download(scriptPos, redownload=False):
    if redownload:
        mkdirs("discarded_db")
        if os.path.exists("%s/db"%scriptPos):
            cmd = "mv %s/db discarded_db"%scriptPos
            obj = Popen(cmd,shell=True)
            obj.wait()

    if not os.path.exists(os.path.join(scriptPos,"db")):
        print("db not exist! download database ...")
        #url="https://zenodo.org/record/6975142/files/db_v0.2.3.tar.gz"
        url="https://zenodo.org/record/8101942/files/db_v0.3.1.tar.gz"
        file=os.path.split(url)[-1]
        cmd = "wget {0} -P {1} && cd {1} && tar -zxvf {2} && rm -rf {2}".format(url, scriptPos, file)
        obj = Popen(cmd,shell=True, stdout=PIPE)
        obj.wait()

        check_db_md5(scriptPos)
    else:
        print("db exist!")


def bayes_classifier_single(inputfile, prefix, wd,
        db_name="prokaryte",
        hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-3,
        blastp_para="-num_threads 3",
        hmmer_para="--noali --cpu 3",
        mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 3"):
    '''
    Aim: single predict lifestyle

    Process: have 2 phase
        phase 1: align to the integrase and excisionase(16 pfam id in total) -- hmmer
        phase 2: align to our protein database -- mmseqs easy-search
        phase 3: dected innovirues -- hmmer and blastp
        if temperate appear in either one phase --> final will be temperate

    Usage: bayes_classifier_single(inputfile,prefix,wd,hmm_creteria=1e-5,mmseqs_creteria=1e-5)
        inputfile: protein set from a bin|a genome.
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)

    Return:[prefix,inte_label,excision_label,pfam_label,p_total_temperate,p_total_lytic,bc_label,final_label]
    '''
    print("## %s start!"%prefix)
    mkdirs(wd)
    #pfam_label, bc_label, final_label = "Lytic", "Lytic", "Lytic"
    pfam_label, bc_label, final_label = "Virulent", "Virulent", "Virulent"
    fileDir = os.path.dirname(os.path.abspath(__file__))


    # phase 1 hmmsearch for pfam
    integrase_hmm = os.path.join(fileDir, "db/integrase_pfv34.hmm")
    excisionase_hmm = os.path.join(fileDir, "db/excisionase_pfv34.hmm")

    # run hmmsearch for integrease and excisionase
    pfam_wd = os.path.join(wd, "BC_pfam")
    mkdirs(pfam_wd)
    inte_opt = runHmmsearch(inputfile, "%s.BC_integrase" %prefix,
                            pfam_wd, integrase_hmm, hmmer_para)
    excision_opt = runHmmsearch(inputfile, "%s.BC_excisionase" % prefix,
                            pfam_wd, excisionase_hmm, hmmer_para)
    # print(inte_opt,excision_opt)

    # load pfam result
    inte_label, excision_label = [0, 0]
    inte_annoD = load_hmmsearch_opt(inte_opt, creteria=hmm_creteria)
    excision_annoD = load_hmmsearch_opt(excision_opt, creteria=hmm_creteria)

    if inte_annoD:
        inte_label = len(inte_annoD)
    if excision_annoD:
        excision_label = len(excision_annoD)

    if inte_label or excision_label:
        #pfam_label = "Lysogenic"
        pfam_label = "Temperate"
    # print(prefix,inte_label,excision_label,pfam_label)

    # phase 2 bayes classifier
    # run mmseqs search
    if db_name == "all":
        print("Using prokaryote and eukaryote protein as DB")
        bc_mmseqsDB = os.path.join(fileDir, "db/bayes_mmseqs_index/all.protein")
    elif db_name == "prokaryote":
        print("Using only prokaryote protein as DB")
        bc_mmseqsDB = os.path.join(fileDir, "db/bayes_mmseqs_index/prokaryote_only")
    else:
        print("Please check parameter -D")

    mmseqs_wd = os.path.join(wd, "BC_mmseqs")
    mmseqs_prefix = "%s.BC_mmseqs" % prefix
    mmseq_opt = runMmseqsEasysearch(inputfile, mmseqs_prefix, mmseqs_wd, bc_mmseqsDB,
                                    otherPara = mmseqs_para)
    # print(mmseq_opt)

    # load score file
    if db_name == "all":
        score_file = os.path.join(fileDir, "db/all_mmseq2_linclust_cluster.stat.scoreOpt.tsv")
    elif db_name == "prokaryote":
        score_file = os.path.join(fileDir, "db/prokaryote_only_mmseq2_linclust_cluster.stat.scoreOpt.tsv")
    else:
        print("Please check parameter -D")

    member2scoreD = load_scoreD(score_file)

    # read mmseqs easy search file
    p_total_temperate, p_total_lytic, bc_label, match_gene_number = calcaulate_score(
        mmseq_opt, member2scoreD, creteria = mmseqs_creteria)
    print(prefix, inte_label, excision_label, pfam_label,
          p_total_temperate, p_total_lytic, bc_label)

    # phase 2B combine result pfam and bayes classifier
    #if pfam_label == "Lysogenic" or bc_label == "Lysogenic":
    #    final_label = "Lysogenic"
    if pfam_label == "Temperate" or bc_label == "Temperate":
        final_label = "Temperate"


    # phase 3: detect Innovirues
    # way1: based on marker gene
    inno_hmmDB = os.path.join(fileDir, "db/Final_marker_morph.hmm")
    inno_blastPre = os.path.join(fileDir, "db/Marker_ALV1")
    inno_dect_wd = os.path.join(wd, "BC_Inno")
    inno_prefix = "%s.Inno" % prefix
    ### in this 1e-3 is used to identify the PI-like
    inno_res = inoviruses_PI_like_gene_search(inputfile, inno_prefix, inno_dect_wd, inno_hmmDB, inno_blastPre,
                              hmm_evalue=blastp_creteria, blastp_evalue=blastp_creteria,
                              blastp_para=blastp_para, hmmer_para=hmmer_para)
    if inno_res == True: final_label = "Chronic"
    print("## %s end!"%prefix)
    return [prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label, match_gene_number]


#def bayes_classifier_batch(inputfile, wd, summaryfile="BC_predict.summary",
#         hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-3,
#         blastp_para="-num_threads 3",
#         hmmer_para="--noali --cpu 3",
#         mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
#                     --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):
#
#    '''
#    Aim: batch predict lifestyle
#    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5,mmseq_creteria=1e-5)
#        inputfile: a tab seperate file contain two column.
#            first column: sample name;
#            second column: path of the protein file;
#        wd: work path where put the result
#        summaryfile: file name for summary of the predict output. location will be under the wd path
#        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
#        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
#        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
#    '''
#
#    mkdirs(wd)
#    opt = open(os.path.join(wd, summaryfile), "w")
#    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
#    opt.write(header)
#
#    with open(inputfile) as f:
#        for line in f:
#            sample_name, path = line.strip("\n").split("\t")
#            res = bayes_classifier_single(
#                path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
#            #prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label = res
#            res.append(path)
#            opt.write("\t".join([str(i) for i in res])+"\n")
#            opt.flush()
#    opt.close()
#
#
#def bayes_classifier_contig(inputfile, wd, summaryfile="BC_predict.summary",
#          hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-3,
#          blastp_para="-num_threads 3",
#          hmmer_para="--noali --cpu 3",
#          mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
#                      --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):
#
#    '''
#    Aim: predict lifestyle for contig file
#    (WARNING: treat each *seq* as an indepent virus sequence, not suit for a genome contain multiple sequence)
#
#    Process: prodigal --> lifestyleBC
#
#    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5,mmseq_creteria=1e-5)
#        inputfile: contig file
#        wd: work path where put the result
#        summaryfile: file name for summary of the predict output. location will be under the wd path
#        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
#        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
#        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
#    '''
#
#    mkdirs(wd)
#    opt = open(os.path.join(wd, summaryfile), "w")
#    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
#    opt.write(header)
#
#    for seq in SeqIO.parse(inputfile,"fasta"):
#        tmpfile="./%s.tmp"%seq.id
#        SeqIO.write(seq,tmpfile,"fasta")
#        faaFile = runProdigal(tmpfile, seq.id, "%s/BC_prodigal"%wd, program="meta", otherPara="-g 11")
#        res = ["sample_name"] + ['NA']*6 + [0]
#        if os.path.getsize(faaFile) != 0:
#            res = bayes_classifier_single(
#                path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
#                #  faaFile, seq.id, wd, hmm_creteria, mmseqs_creteria)
#        #prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label, match_gene_number = res
#        res.append(faaFile)
#        opt.write("\t".join([str(i) for i in res])+"\n")
#        opt.flush()
#        os.remove(tmpfile)
#    opt.close()
#
#
#def bayes_classifier_genomes(inputfile, wd, summaryfile="BC_predict.summary",
#           hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-3,
#           blastp_para="-num_threads 3",
#           hmmer_para="--noali --cpu 3",
#           mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
#                       --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):
#    '''
#    Aim: predict lifestyle for one genome file which contain multiple seq
#    (WARNING: treat *inputfile* as an indepent virus sequence(no matter how many seq are there), suit for a genome contain multiple sequence)
#
#    Process: prodigal --> lifestyleBC
#
#    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-5)
#        inputfile: a tab seperate file contain two column.
#            first column: sample name;
#            second column: path of the genome file;
#        wd: work path where put the result
#        summaryfile: file name for summary of the predict output. location will be under the wd path
#        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
#        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
#        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
#    '''
#
#    mkdirs(wd)
#    opt = open(os.path.join(wd, summaryfile), "w")
#    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
#    opt.write(header)
#
#    with open(inputfile) as f:
#        for line in f:
#            sample_name, path = line.strip("\n").split("\t")
#            faaFile = runProdigal(path, sample_name, "%s/BC_prodigal"%wd, program="meta", otherPara="-g 11")
#            res = [sample_name] + ['NA']*6 + [0]
#            if os.path.getsize(faaFile) != 0:
#                res = bayes_classifier_single(
#                     # faaFile, sample_name, wd, pfam_creteria, mmseqs_creteria)
#                    path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
#            #prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label, match_gene_number = res
#            res.append(faaFile)
#            opt.write("\t".join([str(i) for i in res])+"\n")
#            opt.flush()
#    opt.close()
#
#
if __name__ == "__main__":
#    # for single predict
#    #bayes_classifier_single("./example/simulate_art_sample1.21.faa", "simulate_art_sample1.21", "./test")
#    # bayes_classifier_single("./example/simulate_art_sample1.1.faa","simulate_art_sample1.1","./test")
#    # bayes_classifier_single("./example/simulate_art_sample1.5.faa","simulate_art_sample1.5","./test")
#    # for batch predict
#    bayes_classifier_batch("./example/example.list","./batch_test","BC_predict.summary")
    scriptPos = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(scriptPos)
