#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor, as_completed
from Replidec.Replidec import bayes_classifier_single, runProdigal, checkdb_and_download
from Replidec.utility import mkdirs
import os
from Bio import Seq,SeqIO
import time
from subprocess import Popen



def bayes_classifier_batch(inputfile, wd, db_name="prokaryte", summaryfile="BC_predict.summary",threads=10,
         hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-5,
         blastp_para="-num_threads 3",
         hmmer_para="--noali --cpu 3",
         mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
                     --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):

    '''
    Aim: batch predict lifestyle
    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5,mmseq_creteria=1e-5, blastp_creteria=1e-5)
        inputfile: a tab seperate file contain two column.
            first column: sample name;
            second column: path of the protein file;
        wd: work path where put the result
        summaryfile: file name for summary of the predict output. location will be under the wd path
        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
    '''
    print("Check db")
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)

    mkdirs(wd)
    opt = open(os.path.join(wd, summaryfile), "w")

    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
    opt.write(header)

    executor = ThreadPoolExecutor(max_workers=threads)
    all_task = []
    kwargsD={"hmm_creteria":hmm_creteria, "mmseqs_creteria":mmseqs_creteria, "blastp_creteria":blastp_creteria,
             "blastp_para":blastp_para,"hmmer_para":hmmer_para,"mmseqs_para":mmseqs_para,
             "db_name":db_name}

    ## n is try to control the sceond job will be submit after the first, cause some env or dir will be confilct when run all jobs at same time
    n=0
    with open(inputfile) as f:
        for line in f:
            n+=1
            sample_name, path = line.strip("\n").split("\t")
            if n==2:
                time.sleep(10)
            all_task.append(executor.submit(bayes_classifier_single, path, sample_name, wd, **kwargsD))
            #res = bayes_classifier_single(
                #path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
            #prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label = res
        for future in as_completed(all_task):
            res = future.result()
            res.append(path)
            opt.write("\t".join([str(i) for i in res])+"\n")
            #opt.flush()
    opt.close()

def bayes_classifier_contig(inputfile, wd, db_name="prokaryte",summaryfile="BC_predict.summary",threads=10,
          hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-5,
          blastp_para="-num_threads 3",
          hmmer_para="--noali --cpu 3",
          mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
                      --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):

    '''
    Aim: predict lifestyle for contig file
    (WARNING: treat each *seq* as an indepent virus sequence, not suit for a genome contain multiple sequence)

    Process: prodigal --> lifestyleBC

    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5,mmseq_creteria=1e-5, blastp_creteria=1e-5)
        inputfile: contig file
        wd: work path where put the result
        summaryfile: file name for summary of the predict output. location will be under the wd path
        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
    '''
    print("Check db")
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)

    mkdirs(wd)
    opt = open(os.path.join(wd, summaryfile), "w")
    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
    opt.write(header)

    executor = ThreadPoolExecutor(max_workers=threads)
    all_task = []
    kwargsD={"hmm_creteria":hmm_creteria, "mmseqs_creteria":mmseqs_creteria, "blastp_creteria":blastp_creteria,
              "blastp_para":blastp_para,"hmmer_para":hmmer_para,"mmseqs_para":mmseqs_para,
              "db_name":db_name}

    faaDict = {}
    for seq in SeqIO.parse(inputfile,"fasta"):
        tmpdir = "%s/tmp"%wd
        mkdirs(tmpdir)
        seq.id=seq.id.replace("|","_")
        tmpfile="%s/%s.tmp"%(tmpdir,seq.id)
        SeqIO.write(seq,tmpfile,"fasta")
        faaFile = runProdigal(tmpfile, seq.id, "%s/BC_prodigal"%wd, program="meta", otherPara="-g 11")
        faaDict[seq.id] = faaFile
    #print(faaDict)

    n = 0
    for sample_name,faaFile in faaDict.items():
        if os.path.getsize(faaFile) != 0:
            n+=1
            if n==2:
                time.sleep(10)
            all_task.append(executor.submit(bayes_classifier_single, faaFile, sample_name, wd, **kwargsD))
            ## 串行(deprecated)
            #res = bayes_classifier_single(
            #    path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
    for future in as_completed(all_task):
        res = future.result()
        faaFile = faaDict[res[0]]
        tmpfile_com="%s/%s.tmp"%(tmpdir,res[0])
        res.append(faaFile)
        opt.write("\t".join([str(i) for i in res])+"\n")
        #opt.flush()
        os.remove(tmpfile_com)
    opt.close()


def bayes_classifier_genomes(inputfile, wd, db_name="prokaryte",summaryfile="BC_predict.summary",threads=10,
           hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-5,
           blastp_para="-num_threads 3",
           hmmer_para="--noali --cpu 3",
           mmseqs_para="-s 7 --max-seqs 1 --alignment-mode 3\
                       --alignment-output-mode 0 --min-aln-len 40 --cov-mode 0 --greedy-best-hits 1 --threads 3"):
    '''
    Aim: predict lifestyle for one genome file which contain multiple seq
    (WARNING: treat *inputfile* as an indepent virus sequence(no matter how many seq are there), suit for a genome contain multiple sequence)

    Process: prodigal --> lifestyleBC

    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,hmm_creteria=1e-5, mmseqs_creteria=1e-5, blastp_creteria=1e-5)
        inputfile: a tab seperate file contain two column.
            first column: sample name;
            second column: path of the genome file;
        wd: work path where put the result
        summaryfile: file name for summary of the predict output. location will be under the wd path
        hmm_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
        blastp_creteria: creteria to filter blastp evalue greater than x (default: 1e-5)
    '''
    print("Check db")
    fileDir = os.path.dirname(os.path.abspath(__file__))
    checkdb_and_download(fileDir)

    mkdirs(wd)
    opt = open(os.path.join(wd, summaryfile), "w")
    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tmatch_gene_number\tpath\n"
    opt.write(header)

    executor = ThreadPoolExecutor(max_workers=threads)
    all_task = []
    kwargsD={"hmm_creteria":hmm_creteria, "mmseqs_creteria":mmseqs_creteria, "blastp_creteria":blastp_creteria,
              "blastp_para":blastp_para,"hmmer_para":hmmer_para,"mmseqs_para":mmseqs_para,
              "db_name":db_name}
    faaDict = {}

    with open(inputfile) as f:
        for line in f:
            sample_name, path = line.strip("\n").split("\t")
            faaFile = runProdigal(path, sample_name, "%s/BC_prodigal"%wd, program="meta", otherPara="-g 11")
            #res = [sample_name] + ['NA']*6 + [0]
            faaDict[sample_name] = faaFile

        n = 0
        for sample_name,faaFile in faaDict.items():
            if os.path.getsize(faaFile) != 0:
                n+=1
                if n==2:
                    time.sleep(10)
                all_task.append(executor.submit(bayes_classifier_single, faaFile, sample_name, wd, **kwargsD))
                ## 串行(deprecated)
                #res = bayes_classifier_single(
                #    path, sample_name, wd, hmm_creteria, mmseqs_creteria,blastp_creteria,blastp_para,hmmer_para,mmseqs_para)
        for future in as_completed(all_task):
            res = future.result()
            faaFile = faaDict[res[0]]
            res.append(faaFile)
            opt.write("\t".join([str(i) for i in res])+"\n")
            #opt.flush()
    opt.close()

if __name__ == "__main__":
    ##tested
    bayes_classifier_batch("./example/example.list","./batch_test",summaryfile="BC_predict.summary")
    #bayes_classifier_contig("./example/test.contig.small.fa", "./contig_test", summaryfile="BC_predict.summary",threads=10)
    #bayes_classifier_genomes("example/genome_test.index", "./genome_test", summaryfile="BC_predict.summary",threads=10)
