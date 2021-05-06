from Bio import SeqIO
import sys

def divide_fasta_into_ref_and_query(ref_query, ref):
    
    with    open(ref,'r')                as ihandle       :
        refseq_itr         = SeqIO.parse(ihandle, "fasta")
        refseq_name_set    = set()
        for record in refseq_itr:
            refseq_name_set.add(record.name)

    with    open(ref_query,'r')          as ihandle       ,\
            open(ref_query+".ref",'w')   as ohandle_ref   ,\
            open(ref_query+".query",'w') as ohandle_query :
        ref_query_itr      = SeqIO.parse(ihandle,"fasta")
        ref_query_name_set = {}
        for record in ref_query_itr:
            if(record.name in refseq_name_set):
                SeqIO.write(
                    record        ,
                    ohandle_ref   ,
                    "fasta"
                    )
            else:
                SeqIO.write(
                    record        ,
                    ohandle_query ,
                    "fasta"
                    )

if __name__ == "__main__":
    #main function
    argvs=sys.argv
    divide_fasta_into_ref_and_query(argvs[1],argvs[2])