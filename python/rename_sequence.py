from Bio import SeqIO
import random
import sys
import shutil
import os
import subprocess

def rename_sequence(in_fname,out_fname):
    with open(in_fname) as origin, open(out_fname, 'w') as renamed:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k = 0
        for s in input_itr:
            s.id = "s"+str(k)
            s.description = "s"+str(k)
            k = k +1
            SeqIO.write(s, renamed, "fasta")


def outgroup_check_fast(in_fname):
    exist_root = False
    with open(in_fname) as origin:
        for line in origin:
            if (line           == ">root\n"): exist_root = True; break
            if (line.split()[0]== "root"   ): exist_root = True; break
    return exist_root

def count_sequence(in_fname):
    with open(in_fname) as origin:
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

def count_sequence_fast(in_fname):
    with open(in_fname) as handle:
        k, l=0, 0
        for line in handle:
            if(line[0]==">"): k+=1
            elif(k==1): l+=len(line)-1
    return [k,l] # k: number of sequence, n: sequence length of first sequence (outgroup)

def random_sampling(in_fname,out_fname,subsample_size,seed,n=None, file_format = "fasta"):
    if(seed=="random"): random.seed(int(random.randint(0,99999)))
    elif(len(seed)!=0):random.seed(int(seed))
    else:print("-r Error: invalid random seed!")
    if (n==None):
        n=count_sequence_fast(in_fname)[0]
    if (n > subsample_size):
        rand_idx=random.sample(range(n-1),subsample_size)
        rand_idx.sort()
    else:
        rand_idx=list(range(n-1))

    if ( file_format == "fasta" ):
        with open(out_fname, 'w') as subs:
            with open(in_fname) as allseq:
                allseq_itr = SeqIO.parse(allseq, "fasta")
                for s in allseq_itr:
                    if(s.id=="root"):
                        SeqIO.write(s, subs, "fasta")
            
            added_seqs = set()
            with open(in_fname) as allseq:
                allseq_itr = SeqIO.parse(allseq, "fasta")
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
                            i += 1
                        k += 1
    elif (file_format=="edit"):
        with open(in_fname, 'r') as rhandle, open(out_fname, 'w') as whandle:
            k = 0
            edits_str_set = set()
            for line in rhandle:
                if ( name != "root"):
                    name      = line.split()[0]
                    if ( k==rand_idx[i] ):
                        edits_str = line.split()[1]
                        if edits_str not in edits_str_set:
                            whandle.write(name + "\t" + edits_str + "\n")
                            edits_str_set.add(edits_str)
                        i += 1
                    k += 1
                

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    rename_sequence(argvs[1],argvs[2])