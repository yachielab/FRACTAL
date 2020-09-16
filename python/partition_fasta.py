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

def classify_sequences(inputFASTA_filehandle, seqname2handle, outhandle_list):
    print(seqname2handle)
    # open input FASTAfile
    records = SeqIO.parse(inputFASTA_filehandle, "fasta")
    
    for record in records:
        print (record.name)
        if (record.name == "root"):
            for outhandle in outhandle_list:
                SeqIO.write(record, outhandle, "fasta")
        else:
            outhandle = seqname2handle[record.id]
            SeqIO.write(record, outhandle, "fasta")

def partition_fasta(inputFASTA_filepath, outputFASTA_dirpathlist, seqname2dirpath):
    is_gzipped = (inputFASTA_filepath.split(".")[-1] == "gz")
    if is_gzipped:
        ist  = gzip.open(inputFASTA_filepath, 'rt')
        ist2 = gzip.open(inputFASTA_filepath,'rt')
    else:
        ist  = open(inputFASTA_filepath, 'r')
        ist2 = open(inputFASTA_filepath, 'r')
    ost_list = []
    dirpath2handle = {}
    for outputFASTA_dirpath in outputFASTA_dirpathlist:
        outputFASTA_filepath = outputFASTA_dirpath + "/" + inputFASTA_filepath.split("/")[-1].split(".gz")[0]
        ost_list.append(open(outputFASTA_filepath, 'w'))
        dirpath2handle[outputFASTA_dirpath] = ost_list[-1]

    seqname_set = set()
    records = SeqIO.parse(ist, "fasta")
    for record in records:
        seqname_set.add(record.name)
 
    seqname2handle = {}
    with open(seqname2dirpath, 'r') as dicst:
        for line in dicst:
            line     = line.split("\n")[0]
            seqname  = line.split("\t")[0]
            dirpath  = line.split("\t")[1]
            if seqname in seqname_set:
                seqname2handle[seqname] = dirpath2handle [ dirpath ]

    classify_sequences(ist2, seqname2handle, ost_list)

    ist.close()
    ist2.close()
    for ost in ost_list:
        ost.close()

if __name__ == "__main__":
    partition_fasta(
        inputFASTA_filepath     = sys.argv[1], 
        outputFASTA_dirpathlist = sys.argv[2].split(":"), 
        seqname2dirpath         = sys.argv[3]
        )