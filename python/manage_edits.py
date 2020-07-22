from Bio import SeqIO, Seq
import sys
import gzip

def edit2editlist(edit_file):
    
    inhandle = gzip.open(edit_file, 'rt')
    
    edit2seqcount  = {}
    total_seqcount = 0
    for line in inhandle:
        total_seqcount+= 1
        seq_name      = line.split()[0]
        if (len(line.split()) > 1):
            seq_edit_list = line.split()[1].split(";")
        else:
            seq_edit_list = []
        for seq_edit in seq_edit_list:
            if seq_edit in edit2seqcount.keys():
                edit2seqcount[seq_edit] += 1
            else:
                edit2seqcount[seq_edit] = 1

    editlist = []
    for key_edit in edit2seqcount:
        if (edit2seqcount[key_edit] < total_seqcount):
            editlist.append(key_edit)
            editlist.sort()

    inhandle.close()

    return editlist

def edit2fasta(edit_file, edit_list, out_gz = True):
    inhandle = gzip.open(edit_file, 'rt')
    outhandle = gzip.open( edit_file + ".fa.gz", 'wt' )
    
    if (True):
        for line in inhandle:
            seq_name      = line.split()[0]
            if len(line.split()) > 1:
                seq_edit_list = line.split()[1].split(";")
            else:
                seq_edit_list = []
            seq_edit_set  = set(seq_edit_list)
            seq_str       = ""
            seq_half_pos  = len(edit_list)//2 + 1
            for edit in edit_list[:seq_half_pos]:
                if (edit in seq_edit_set):
                    seq_str = seq_str + "T"
                else:
                    seq_str = seq_str + "C"
            for edit in edit_list[seq_half_pos:]:
                if (edit in seq_edit_set):
                    seq_str = seq_str + "A"
                else:
                    seq_str = seq_str + "G"
            outhandle.write(">"+seq_name+"\n"+seq_str+"\n")
            
    inhandle.close()
    outhandle.close()

'''
command line argument: "<editsets .edit file path> <edit list .txt file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    
    editlist_file = argvs[2]
    editlist = []
    with open(editlist_file, 'r') as handle:
        for line in handle:
            editlist.append(line.split("\n")[0])

    edit2fasta(argvs[1],editlist,out_gz=True)