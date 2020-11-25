from Bio import SeqIO
import random
import sys
import shutil
import os
import subprocess
import gzip

def almighty_open(fpath, mode='r'):
    if fpath.split(".")[-1] == "gz":
        return gzip.open(fpath, mode+'t')
    else:
        return open(fpath, mode)

def rename_sequence(in_fname,out_fname,file_format='fa'):
    is_gzipped = (in_fname.split(".")[-1] == "gz")
    if (is_gzipped):
        renamed   = gzip.open(out_fname, 'wt')
        origin    = gzip.open(in_fname, 'rt')
    else:
        renamed   = open(out_fname, 'w')
        origin    = open(in_fname, 'r')

    if   file_format == 'fa':
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k = 0
        name2renamed = {"root":"s0"}
        for s in input_itr:
            name2renamed[s.name] = "s"+str(k)
            s.id = "s"+str(k)
            s.description = "s"+str(k)
            SeqIO.write(s, renamed, "fasta")
            k += 1
    elif file_format == 'edit':
        name2renamed = {"root":"s0"}
        renamed.write("s0\t\n")
        k = 1
        for line in origin:
            if line.split()[0]!="root":
                name     = line.split()[0]
                if (len(line.split())==2): edit_str = line.split()[1]
                name2renamed[name] = "s"+str(k)
                if (len(line.split())==2): renamed.write("s"+str(k)+"\t"+edit_str+"\n")
                elif (len(line.split())==1): renamed.write("s"+str(k)+"\t\n")
                k += 1

    renamed.close()
    origin.close()

    return name2renamed

def outgroup_check_fast(in_fpathlist, file_format):
    exist_root = False
    for in_fname in in_fpathlist:
        if (exist_root): break
        is_gzipped = (in_fname.split(".")[-1] == "gz")
        
        # file open
        if ( is_gzipped ):
            origin = gzip.open(in_fname, 'rt') 
        else:
            origin = open(in_fname, 'r') 

        idx = 0
        if (file_format == "fasta"):
            records = SeqIO.parse(origin, file_format)
            for record in records:
                if (record.name  == "root"):
                    exist_root = True
                    root_fpath = "/".join(in_fname.split("/")[:-1]) + "/root.fa"
                    SeqIO.write(record, root_fpath, 'fasta')
                    break
        elif (file_format == "edit"):
            for line in origin:
                if   (line.split()[0]== "root"): 
                    exist_root = True
                    root_fpath = in_fname.split("/")[-1] + "/root.edit"
                    with open(root_fpath, 'w') as root_handle:
                        root_handle.write(line)
                    break
                else: 
                    idx += 1

    # file close
    origin.close()
    if (exist_root):
        return root_fpath
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
def count_sequence_fast(in_fpathlist, form = "fa"): # form: fa or edit
    fpath2seqcount = {}
    for in_fname in in_fpathlist:
        is_gzipped = (in_fname.split(".")[-1] == "gz")

        if (is_gzipped):    gunzip = "| gunzip"
        else:               gunzip = ""

        if (form == "fa"):
            grep_command = "| grep '>'"
        else:
            grep_command = ""
        seq_count_str = (
            subprocess.Popen(
                "cat " + in_fname + " " + gunzip + grep_command + "| wc -l",
                stdout=subprocess.PIPE,
                shell=True
                ).communicate()[0]
                ).decode('utf-8')
        fpath2seqcount[in_fname] = int(seq_count_str)
    return fpath2seqcount

def random_sampling( # fasta only
    in_dirpath,
    out_fname,
    subsample_size,
    fpath2seqcount,
    root_fpath,
    total_seqcount=None, 
    file_format   ="fa",
    in_fpath      =None
    ):

    # set input file path list
    if (in_fpath == None):
        fpath_list     = sorted([in_dirpath + "/" + fname for fname in os.listdir(in_dirpath)])
    else:
        fpath_list     = [in_fpath]

    
    # count total number of sequences
    if (total_seqcount == None):
        fpath2seqcount = count_sequence_fast(fpath_list,form=file_format)
        total_seqcount = 0
        for fpath in fpath_list:
            seqcount              = fpath2seqcount[fpath]
            total_seqcount       += seqcount

    # get sequence indices to extract
    if (total_seqcount > subsample_size):
        rand_idx_list=random.sample(range(total_seqcount),subsample_size) # n-1: not include root
        rand_idx_list.sort()
    else:
        rand_idx_list=list(range(total_seqcount))
    rand_idx_set = set(rand_idx_list)

    # get sequence index > file name + local index
    idx2fpath_localindex = {}
    idx            = 0
    fpath2localidx = {}
    for fpath in fpath_list:
        for local_idx in range(fpath2seqcount[fpath]):
            if (idx in rand_idx_set):
                try   : fpath2localidx[fpath].append(local_idx)
                except: fpath2localidx[fpath]     = [local_idx]
            idx += 1
    
    # sequence extraction
    seq_set      = set()
    seqname_list = []
    with open(out_fname, 'w') as ost:
        print(file_format)
        if (file_format == 'fa'):
            print(root_fpath)
            # write root sequence
            with open(root_fpath, 'r') as ist:
                records        = SeqIO.parse(ist, 'fasta')
                for record in records:
                    seq_set.add(str(record.seq))
                    SeqIO.write(record, ost, 'fasta')
                    seqname_list.append(record.name)
            for fpath in sorted(list(fpath2localidx.keys())):
                with almighty_open(fpath, 'r') as ist:
                    local_idx_list = fpath2localidx[fpath]
                    records        = SeqIO.parse(ist, 'fasta')
                    i = 0
                    for record in records:
                        if len(local_idx_list) < 1:
                            break 
                        if i == local_idx_list[0]:
                            if ( str(record.seq) not in seq_set ):
                                seq_set.add(str(record.seq))
                                SeqIO.write(record, ost, 'fasta')
                                seqname_list.append(record.name)
                            local_idx_list.pop(0)
                        i += 1
        if (file_format == 'edit'):
            ost.write("root\t\n")
            for fpath in sorted(list(fpath2localidx.keys())):
                with almighty_open(fpath, 'r') as ist:
                    local_idx_list = fpath2localidx[fpath]
                    i = 0
                    for line in ist:
                        if len(local_idx_list) < 1:
                            break 
                        if i == local_idx_list[0]:
                            seq_set.add(line)
                            ost.write(line)
                            seqname_list.append(line.split()[0])
                            local_idx_list.pop(0)
                        i += 1
    return seqname_list

# To be deleted
def random_sampling_fasta( # fasta only
    in_dirpath,
    out_fname,
    subsample_size,
    fpath2seqcount,
    root_fpath,
    total_seqcount=None, 
    file_format   ="fa",
    in_fpath      =None
    ):

    # set input file path list
    if (in_fpath == None):
        fpath_list     = sorted([in_dirpath + "/" + fname for fname in os.listdir(in_dirpath)])
    else:
        fpath_list     = [in_fpath]
    
    # count total number of sequences
    if (total_seqcount == None):
        fpath2seqcount = count_sequence_fast(fpath_list)
        total_seqcount = 0
        for fpath in fpath_list:
            seqcount              = fpath2seqcount[fpath]
            total_seqcount       += seqcount

    # get sequence indices to extract
    if (total_seqcount > subsample_size):
        rand_idx_list=random.sample(range(total_seqcount),subsample_size) # n-1: not include root
        rand_idx_list.sort()
    else:
        rand_idx_list=list(range(total_seqcount))
    rand_idx_set = set(rand_idx_list)

    # get sequence index > file name + local index
    idx2fpath_localindex = {}
    idx            = 0
    fpath2localidx = {}
    for fpath in fpath_list:
        for local_idx in range(fpath2seqcount[fpath]):
            if (idx in rand_idx_set):
                try   : fpath2localidx[fpath].append(local_idx)
                except: fpath2localidx[fpath]     = [local_idx]
            idx += 1
    
    # sequence extraction
    seq_set      = set()
    seqname_list = []
    with open(out_fname, 'w') as ost:
        if (file_format == 'fa'):
            # write root sequence
            with open(root_fpath, 'r') as ist:
                records        = SeqIO.parse(ist, 'fasta')
                for record in records:
                    seq_set.add(str(record.seq))
                    SeqIO.write(record, ost, 'fasta')
                    seqname_list.append(record.name)
            for fpath in sorted(list(fpath2localidx.keys())):
                with almighty_open(fpath, 'r') as ist:
                    local_idx_list = fpath2localidx[fpath]
                    records        = SeqIO.parse(ist, 'fasta')
                    i = 0
                    for record in records:
                        if len(local_idx_list) < 1:
                            break 
                        if i == local_idx_list[0]:
                            if ( str(record.seq) not in seq_set ):
                                seq_set.add(str(record.seq))
                                SeqIO.write(record, ost, 'fasta')
                                seqname_list.append(record.name)
                            local_idx_list.pop(0)
                        i += 1
        if (file_format == 'edit'):
            for fpath in sorted(list(fpath2localidx.keys())):
                with almighty_open(fpath, 'r') as ist:
                    local_idx_list = fpath2localidx[fpath]
                    i = 0
                    for line in ist:
                        if len(local_idx_list) < 1:
                            break 
                        if i == local_idx_list[0]:
                            seq_set.add(line)
                            ost.write(line, ost)
                            seqname_list.append(line.split()[0])
                            local_idx_list.pop(0)
                        i += 1
    return seqname_list

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    rename_sequence(argvs[1],argvs[2])