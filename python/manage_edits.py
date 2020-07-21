from Bio import SeqIO


def edit2editlist(edit_file):
    
    is_gzipped = (edit_file.split(".")[-1] == "gz")
    if ( is_gzipped ):
        inhandle = gzip.open(edit_file, 'rt')
    else:
        inhandle = open(edit_file, 'r')
    
    edit2seqcount  = {}
    total_seqcount = 0
    for line in inhandle:
        total_seqcount+= 1
        seq_name      = line.split()[0]
        seq_edit_list = line.split()[1].split(";")
        for seq_edit in seq_edit_list:
            if seq_edit in edit2seqcount.keys():
                edit2seqcount[seq_edit] += 1
            else:
                edit2seqcount[seq_edit] = 1

    editlist = []
    for key_edit in edit2seqcount:
        if (edit2seqcount[key_edit] == total_seqcount):
            editlist.append(key_edit)
            editlist.sort()

    inhandle.close()

    return editlist

def edit2fasta(edit_file, edit_list):
    is_gzipped = (edit_file.split(".")[-1] == "gz")
    if ( is_gzipped ):
        inhandle = gzip.open(edit_file, 'rt')
    else:
        inhandle = open(edit_file, 'r')

    with open( edit_file + ".fa", 'w' ) as outhandle:

        for line in inhandle:
            seq_name      = line.split()[0]
            seq_edit_list = line.split()[1].split(";")
            seq_edit_set  = set(seq_edit_list)
            seq_str       = ""
            seq_half_pos  = len(edit_list)//2 + 1
            for edit in edit_list[:seq_half_pos]:
                if (edit in edit_list):
                    seq_str = seq_str + "T"
                else:
                    seq_str = seq_str + "C"
            for edit in edit_list[seq_half_pos:]:
                if (edit in edit_list):
                    seq_str = seq_str + "A"
                else:
                    seq_str = seq_str + "G"
            SeqIO.write(seq_str, outhandle, "fasta")
            
    inhandle.close()

'''
command line argument: "<editsets .edit file path> <edit list .txt file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    
    editlist_file = argvs[2]
    editlist = []
    with open(editlist_file, 'r') as handle:
        for line in handle:
            editlist.append(line.split()[0])

    rename_sequence(argvs[1],editlist)