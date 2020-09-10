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

def classify_sequences(inputFASTA_filehandle, seqname2handle):
    # open input FASTAfile
    records = SeqIO.parse(inputFASTA_filehandle, "fasta")
    
    for record in records:
        outhandle = seqname2handle[record.id]
        SeqIO.write(outhandle, record, "fasta")

def partition_fasta(inputFASTA_path, outputFASTA_dirpathlist, seqname2dirpath):

    seqname_set = {}
    with gzip.open(inputFASTA_path, 'rt') as ist:
        records = SeqIO.parse(ist)
        for record in records:
            seqname_set.add(record.name)

    ist = gzip.open(inputFASTA_path,'rt')
    ost_list = []
    dirpath2handle = {}
    for outputFASTA_dirpath in outputFASTA_dirpathlist:
        outputFASTA_filepath = outputFASTA_dirpath + "/" + inputFASTA_path.split("/")[-1]
        ost_list.append(gzip.open(outputFASTA_filepath, 'wt'))
        dirpath2handle[outputFASTA_dirpath] = ost_list[-1]
    
    seqname2handle = {}
    with open(seqname2dirpath, 'r') as dicst:
        for line in dicst:
            line     = line.split("\n")[0]
            seqname  = line.split("\t")[0]
            dirpath = line.split("\t")[1]
            if seqname in seqname_set:
                seqname2handle[seqname] = dirpath2handle [ dirpath ]

    classify_sequences(ist, seqname2handle)

if __name__ == "__main__":
    partition_fasta(
        inputFASTA_path      = sys.argv[1], 
        outputFASTA_pathlist = sys.argv[2].split(":"), 
        seqname2dirpath      = sys.argv[3]
        )