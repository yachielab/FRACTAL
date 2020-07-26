from Bio import SeqIO
import random
import sys
import shutil
import os
import subprocess
import gzip

def rename_sequence(in_fname,out_fname):
    with gzip.open(in_fname, 'rt') as origin, gzip.open(out_fname, 'wt') as renamed:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k = 0
        name2renamed = {}
        for s in input_itr:
            name2renamed[s.name] = "s"+str(k)
            s.id = "s"+str(k)
            s.description = "s"+str(k)
            k = k + 1
            SeqIO.write(s, renamed, "fasta")
    return name2renamed

def outgroup_check_fast(in_fname):
    exist_root = False
    with gzip.open(in_fname, 'rt') as origin:
        for line in origin:
            if (line           == ">root\n"): exist_root = True; break
            if (line.split()[0]== "root"   ): exist_root = True; break
    return exist_root

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

def count_sequence_fast(in_fname):
    with gzip.open(in_fname, 'rt') as handle:
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
        if ( file_format == "fasta" ):
            n=count_sequence_fast(in_fname)[0]
        if ( file_format == "edit" ):
            n = int(subprocess.check_output(['wc', '-l', in_fname]).decode().split(' ')[0])
    if (n > subsample_size):
        rand_idx=random.sample(range(n-1),subsample_size)
        rand_idx.sort()
    else:
        rand_idx=list(range(n-1))

    sample_name_list = ["root"]
    if ( file_format == "fasta" ):
        with gzip.open(out_fname, 'wt') as subs:
            with gzip.open(in_fname, 'rt') as allseq:
                allseq_itr = SeqIO.parse(allseq, "fasta")
                for s in allseq_itr:
                    if(s.id=="root"):
                        SeqIO.write(s, subs, "fasta")    
            
            added_seqs = set()
            with gzip.open(in_fname, 'rt') as allseq:
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
                                sample_name_list.append(s.name)
                            i += 1
                        k += 1
    elif (file_format=="edit"):
        with gzip.open(in_fname, 'rt') as rhandle, gzip.open(out_fname, 'wt') as whandle:
            i = 0
            k = 0
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