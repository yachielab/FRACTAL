from Bio import SeqIO
import random
import sys
import shutil
import os
import subprocess
import gzip

def rename_sequence(in_fname,out_fname):
    is_gzipped = (in_fname.split(".")[-1] == "gz")
    if (is_gzipped):
        renamed   = gzip.open(out_fname, 'wt')
        origin    = gzip.open(in_fname, 'rt')
    else:
        renamed   = open(out_fname, 'w')
        origin    = open(in_fname, 'r')

    input_itr = SeqIO.parse(origin, "fasta")
    # Build a list sequences:
    k = 0
    name2renamed = {"root":"s0"}
    for s in input_itr:
        name2renamed[s.name] = "s"+str(k)
        s.id = "s"+str(k)
        s.description = "s"+str(k)
        k = k + 1
        SeqIO.write(s, renamed, "fasta")

    renamed.close()
    origin.close()

    return name2renamed

def outgroup_check_fast(in_fname, file_format):
    is_gzipped = (in_fname.split(".")[-1] == "gz")

    exist_root = False
    
    # file open
    if ( is_gzipped ):
        origin = gzip.open(in_fname, 'rt') 
    else:
        origin = open(in_fname, 'r') 

    idx = 0
    for line in origin:
        if (file_format == "fasta"):
            if   (line           == ">root\n"): exist_root = True; break
            elif (line[0]== '>'): 
                idx += 1
        if (file_format == "edit"):
            if   (line.split()[0]== "root"   ): exist_root = True; break
            else: 
                idx += 1

    # file close
    origin.close()
    if (exist_root):
        return idx
    else:
        raise Exception('No root')

def count_sequence(in_fname):
    with gzip.open(in_fname, 'rt') as origin:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k=0
        for s in input_itr:
            k+=1
        input_itr = SeqIO.parse(origin, "fasta")
        n=0
        for s in input_itr:
            n+=len(str(s.seq))//k
    return [k,n] # k: number of sequence, n: average sequence length
'''
def count_sequence_fast(in_fname):
    with gzip.open(in_fname, 'rt') as handle:
        k, l=0, 0
        for line in handle:
            if(line[0]==">"): k+=1
            elif(k==1): l+=len(line)-1
    return [k,l] # k: number of sequence, n: sequence length of first sequence (outgroup)
'''
def count_sequence_fast(in_fname):
    is_gzipped = (in_fname.split(".")[-1] == "gz")

    if (is_gzipped):    gunzip = "| gunzip"
    else:               gunzip = ""

    seq_count_str = (
        subprocess.Popen(
            "cat " + in_fname + gunzip + " | grep '>' | wc -l",
            stdout=subprocess.PIPE,
            shell=True
            ).communicate()[0]
            ).decode('utf-8')
    return int(seq_count_str)

def random_sampling_from_splitted( # fasta only
    in_dirname,
    out_fname,
    subsample_size,
    seed,
    Nseq_per_file,
    root_idx,
    n=None, 
    file_format = "fasta",
    is_gzipped = False
    ):

    if (n > subsample_size):
        rand_idx=random.sample(range(1,n),subsample_size) # n-1: not include root
        rand_idx.sort()
    else:
        rand_idx=list(range(1,n))
    
    rand_idx = [root_idx]+rand_idx

    if (is_gzipped):
        gzip_command = "gzip"
    else:
        gzip_command = "cat"

    for k, seq_idx in enumerate(rand_idx):
        file_idx            = seq_idx // Nseq_per_file + 1
        seq_idx_in_the_file = seq_idx %  Nseq_per_file + 1

        file_idx_str        = str(file_idx).zfill(3)
        
        in_fname            = in_dirname+"/"+".".join(in_dirname.split(".")[:-2]).split("/")[-1]+".part_" + file_idx_str + "." + ".".join(in_dirname.split(".")[-2:-1])
        
        command = "seqkit range -r " + str(seq_idx_in_the_file) + ":" + str(seq_idx_in_the_file) + " " + in_fname + "| "+gzip_command+" > " + out_fname+"."+str(k)
        subprocess.call(command, shell = True)
    subprocess.call("cat "+out_fname+".* > " + out_fname+"; rm "+out_fname+".*", shell = True)

    seq_names_str = (
        subprocess.Popen(
            "cat "+out_fname+" | gunzip | grep '>' | tr -d '>' | tr '\n' ','",
            stdout=subprocess.PIPE,
            shell=True
            ).communicate()[0]
            ).decode('utf-8')
    seqnames = seq_names_str.split(",")
    return seqnames
    

def random_sampling(in_fname,out_fname,subsample_size,seed,n=None, file_format = "fasta"):
    is_gzipped = (in_fname.split(".")[-1] == "gz")
    
    if (n==None):
        if ( file_format == "fasta" ):
            n=count_sequence_fast(in_fname)
        if ( file_format == "edit" ):
            n = int((subprocess.Popen('less '+in_fname+' | wc -l', stdout=subprocess.PIPE, shell=True).communicate()[0]).decode('utf-8'))
    if (n > subsample_size):
        rand_idx=random.sample(range(n-1),subsample_size) # n-1: not include root
        rand_idx.sort()
    else:
        rand_idx=list(range(n-1))
    sample_name_list = []
    if ( file_format == "fasta" ):
        # open files
        if (is_gzipped):
            subs   = gzip.open(out_fname, 'wt')
            allseq = gzip.open(in_fname, 'rt')
            allseq2 = gzip.open(in_fname, 'rt')
        else:
            subs   = open(out_fname, 'w')
            allseq = open(in_fname, 'r')
            allseq2 = open(in_fname, 'r')


        allseq_itr = SeqIO.parse(allseq, "fasta")
        for s in allseq_itr:
            if(s.id=="root"):
                SeqIO.write(s, subs, "fasta")    
        added_seqs = set()
        
        allseq_itr = SeqIO.parse(allseq2, "fasta")
        i=0 # index on rand_idx
        k=0 # index on record
        for s in allseq_itr:
            if(i>=len(rand_idx)):
                break
            if(s.id!="root"):
                if(k==rand_idx[i]):
                    if (str(s.seq) not in added_seqs):
                        SeqIO.write(s, subs, "fasta")
                        added_seqs.add(str(s.seq))
                        sample_name_list.append(s.name)
                    i += 1
                k += 1
        # close files
        subs.close()
        allseq.close()
        allseq2.close()
    elif (file_format=="edit"):
        with gzip.open(in_fname, 'rt') as rhandle, gzip.open(out_fname, 'wt') as whandle:
            i = 0
            k = 0 # line number
            edits_str_set = set()
            whandle.write("root\n")
            for line in rhandle:
                name      = line.split()[0]
                if ( name != "root"):
                    if (i>=len(rand_idx)):
                        break
                    if ( k==rand_idx[i] ):
                        if ( len(line.split()) > 1 ):
                            edits_str = line.split()[1]
                        else:
                            edits_str = ""
                        if edits_str not in edits_str_set:
                            whandle.write(name + "\t" + edits_str + "\n")
                            edits_str_set.add(edits_str)
                            sample_name_list.append(name)
                        i += 1
                    k += 1
    return sample_name_list

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    rename_sequence(argvs[1],argvs[2])