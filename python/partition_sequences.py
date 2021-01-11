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

def classify_sequences(inputFASTA_filehandle, seqname2dirpath, dirpath2filepath, is_gzipped, file_format):
    # open input FASTAfile
    outfilepath2Nseq = {}
    filepath2handle  = {}

    if   (file_format=='fa'):   
        records = SeqIO.parse(inputFASTA_filehandle, "fasta")
        
        for record in records:
            if (record.name == "root"):
                None
            else:
                outfilepath = dirpath2filepath[seqname2dirpath[record.id]]

                if (outfilepath not in filepath2handle.keys()):
                    filepath2handle[outfilepath] = open(outfilepath, 'w')
                    outfilepath2Nseq[outfilepath]= 0

                outfilepath2Nseq[outfilepath] += 1
                outhandle   = filepath2handle[outfilepath]
                SeqIO.write(record, outhandle, "fasta")
    elif (file_format=='edit'): 
        for line in inputFASTA_filehandle:
            name         = line.split()[0]
            if len(line.split())==2: edit_pattern = line.split()[1]
            if (name == "root"):
                None
            else:
                outfilepath = dirpath2filepath[seqname2dirpath[name]]
                if (outfilepath not in filepath2handle.keys()):
                    filepath2handle[outfilepath] = open(outfilepath, 'w')
                    outfilepath2Nseq[outfilepath]= 0
                outfilepath2Nseq[outfilepath] += 1
                outhandle = filepath2handle[outfilepath]
                outhandle.write(name + '\t')
                if len(line.split())==2: outhandle.write(edit_pattern + '\n')

    
    for outfilepath in filepath2handle.keys():
        filepath2handle[outfilepath].close()
        if (is_gzipped):
            subprocess.call(
                "gzip " + outfilepath,
                shell=True
                )

        countfilepath = "/".join(outfilepath.split("/")[:-2]) + "/count/" + outfilepath.split("/")[-1] + ".count"
        with open(countfilepath, 'w') as numhandle:
            if (is_gzipped):
                numhandle.write(outfilepath + ".gz\t" + str(outfilepath2Nseq[outfilepath]) + "\n")
            else:
                numhandle.write(outfilepath + "\t"    + str(outfilepath2Nseq[outfilepath]) + "\n")


def partition_sequences(inputFASTA_filepathlist, outputFASTA_dirpathlist, seqname2dir_filepath, file_format = 'fa'):

    data_type = outputFASTA_dirpathlist[0].split("/")[-1]

    seqname_set = set()
    for inputFASTA_filepath in inputFASTA_filepathlist:
        if (os.path.exists(inputFASTA_filepath)):
            is_gzipped = (inputFASTA_filepath.split(".")[-1] == "gz")
            if is_gzipped:
                ist  = gzip.open(inputFASTA_filepath, 'rt')
            else:
                ist  = open(inputFASTA_filepath, 'r')
            if   file_format == 'fa':
                records = SeqIO.parse(ist, "fasta")
                for record in records:
                    seqname_set.add(record.name)
            elif file_format == 'edit':
                for line in ist:
                    seqname_set.add(line.split()[0])
            ist.close()
    if len(seqname_set)==0:
        return

    seqname2dirpath = {}
    with open(seqname2dir_filepath, 'r') as dicst:
        for line in dicst:
            line     = line.split("\n")[0]
            seqname  = line.split("\t")[0]
            dirpath  = line.split("\t")[1]
            if seqname in seqname_set:
                seqname2dirpath[seqname] = dirpath + "/INPUT/" + data_type

    for inputFASTA_filepath in inputFASTA_filepathlist:
        if (os.path.exists(inputFASTA_filepath)):
            is_gzipped = (inputFASTA_filepath.split(".")[-1] == "gz")
            if is_gzipped:
                ist  = gzip.open(inputFASTA_filepath, 'rt')
            else:
                ist  = open(inputFASTA_filepath, 'r')
            filepath2handle = {}
            dirpath2filepath = {}
            for outputFASTA_dirpath in outputFASTA_dirpathlist:
                if inputFASTA_filepath.split(".")[-1] == 'gz':
                    outputFASTA_filepath = outputFASTA_dirpath + "/" + ".".join(inputFASTA_filepath.split("/")[-1].split(".")[:-1])
                else:
                    outputFASTA_filepath = outputFASTA_dirpath + "/" + inputFASTA_filepath.split("/")[-1]
                dirpath2filepath[outputFASTA_dirpath] = outputFASTA_filepath

            classify_sequences(ist, seqname2dirpath, dirpath2filepath, is_gzipped, file_format = file_format)
            ist.close()
            os.remove(inputFASTA_filepath)

if __name__ == "__main__":
    partition_sequences(
        inputFASTA_filepathlist = sys.argv[1].split(":"), 
        outputFASTA_dirpathlist = sys.argv[2].split(":"), 
        seqname2dir_filepath    = sys.argv[3]
        )