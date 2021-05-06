from Bio import SeqIO
import sys


def transform(infile, inform, outform):
    with open(infile,"r") as input_handle, open(infile+"."+outform,"w") as output_handle:
        sequences = SeqIO.parse(input_handle, inform)
        count = SeqIO.write(sequences, output_handle, outform)

'''
main function
        commandline argument: "python3 transform.py <input file path>.fa fasta phylip"
'''
if __name__ == "__main__":
    argvs=sys.argv
    transform(argvs[1], argvs[2], argvs[3])