from Bio import Phylo
from io import StringIO
import sys
import os
import json
import jplace_parse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import random
import gzip

def classify_sequences(inputFASTA_filehandle, outputFASTA_filehandlelist, seqname2dirpath, dirpath2handle):
    # open input FASTAfile
    records = SeqIO.parse(inputFASTA_filehandle, "fasta")
    
    for record in records:
        if (record.name == "root"):
            for outhandle in outputFASTA_filehandlelist:
                SeqIO.write(record, outhandle, "fasta")
        else:
            outhandle = dirpath2handle[seqname2dirpath[record.id]]
            SeqIO.write(record, outhandle, "fasta")

def partition_sequences(inputFASTA_filepathlist, outputFASTA_dirpathlist, seqname2dir_filepath):

    seqname_set = set()
    for inputFASTA_filepath in inputFASTA_filepathlist:
        is_gzipped = (inputFASTA_filepath.split(".")[-1] == "gz")
        if is_gzipped:
            ist  = gzip.open(inputFASTA_filepath, 'rt')
        else:
            ist  = open(inputFASTA_filepath, 'r')
        records = SeqIO.parse(ist, "fasta")
        for record in records:
            seqname_set.add(record.name)
        ist.close()
    seqname2dirpath = {}
    with open(seqname2dir_filepath, 'r') as dicst:
        for line in dicst:
            line     = line.split("\n")[0]
            seqname  = line.split("\t")[0]
            dirpath  = line.split("\t")[1]
            if seqname in seqname_set:
                seqname2dirpath[seqname] = dirpath

    for inputFASTA_filepath in inputFASTA_filepathlist:
        is_gzipped = (inputFASTA_filepath.split(".")[-1] == "gz")
        if is_gzipped:
            ist  = gzip.open(inputFASTA_filepath, 'rt')
        else:
            ist  = open(inputFASTA_filepath, 'r')
        ost_list = []
        dirpath2handle = {}
        for outputFASTA_dirpath in outputFASTA_dirpathlist:
            outputFASTA_filepath = outputFASTA_dirpath + "/" + inputFASTA_filepath.split("/")[-1].split(".gz")[0]
            ost_list.append(open(outputFASTA_filepath, 'w'))
            dirpath2handle[outputFASTA_dirpath] = ost_list[-1]

        classify_sequences(ist, ost_list, seqname2dirpath, dirpath2handle)

        ist.close()
        for ost in ost_list:
            ost.close()
        if (is_gzipped):    
            for outputFASTA_dirpath in outputFASTA_dirpathlist:
                outputFASTA_filepath = outputFASTA_dirpath + "/" + inputFASTA_filepath.split("/")[-1].split(".gz")[0]
                subprocess.call(
                    "gzip " + outputFASTA_filepath,
                    shell=True
                )
        os.remove(inputFASTA_filepath)

if __name__ == "__main__":
    partition_sequences(
        inputFASTA_filepathlist = sys.argv[1].split(":"), 
        outputFASTA_dirpathlist = sys.argv[2].split(":"), 
        seqname2dir_filepath    = sys.argv[3]
        )