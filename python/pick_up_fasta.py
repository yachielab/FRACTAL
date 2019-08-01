import sys
from Bio import SeqIO

# get small fasta file from i*x th sequence to i*(x+1)-1 th sequence of [in_file].fa
def pick_up_fasta(in_file, i, x):
    with open(in_file,'r') as ihandle, open(in_file+"."+str(i), 'w') as ohandle:
        allseq_itr = SeqIO.parse(ihandle, "fasta")
        k=0
        for seq in allseq_itr:
            if(i*x<=k and k<(i+1)*x):
                SeqIO.write(seq,ohandle,'fasta')
            k+=1

if __name__ == "__main__":
    argvs = sys.argv
    print([argvs[1],argvs[2],argvs[3]])
    pick_up_fasta(argvs[1],int(argvs[2]),int(argvs[3]))