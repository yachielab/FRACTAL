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

def partition_fasta(in_fasta_list,num_file,OUT_DIR,wd,jpart,info,treefile,subsamplefa,ROOTING,file_format="fasta"):
    # open .jpart file
    with open(jpart,"r") as jf:
        jp = jf.read()
    # parse json format
    js = json.loads(jp)
    # get tree information
    treestr = js['tree']
    tree = Phylo.read(StringIO(treestr),"newick")
    num_mono = len(tree.get_terminals())

    # get number (decides dir name) and make dir
    num=0
    while(True):
        with open(num_file,'r') as nf:
            try:
                num=int(nf.read())
            except Exception as e:
                print(e)
        with open(num_file,'w') as nf:
            try:
                nf.write(str(num+num_mono))
            except Exception as e:
                print(e)
        try:
            for i in range(num_mono):
                os.mkdir(OUT_DIR+"/d"+str(num+i))
            break
        except Exception as e:
            print(e)
    # create dictionary of leaf -> OUT_DIR/dir path
    DIRdict = {}
    NUMdict = {}
    i = 0
    for leaf in tree.get_terminals():
        DIRdict[leaf.name] = [OUT_DIR+"/d"+str(num+i),0]
        NUMdict[leaf.name.strip('{').strip('}')] = i
        i=i+1

    for fasta_count, in_fasta in enumerate(in_fasta_list):
        ost=[]
        for i in range(num_mono):
            ost.append(gzip.open(OUT_DIR+"/d"+str(num+i)+"/"+in_fasta.split("/")[-1],'wt'))
        para=gzip.open(wd+"/"+in_fasta.split("/")[-1]+".problematic.gz",'wt')

        if (file_format=="fasta"):
            with gzip.open(in_fasta,'rt') as in_handle:
                record = SeqIO.parse(in_handle, "fasta")
                i=0
                for s in record:
                    if(s.id=="root"):
                        for st in ost: #ROOTING=="Origin"
                            SeqIO.write(s, st, "fasta")
                    elif(js["partition"][s.id]=="paraphyletic"):
                        SeqIO.write(s, para, "fasta")
                    else:
                        l = NUMdict[str(js["partition"][s.id])]
                        if( fasta_count == 0 ):
                            DIRdict['{'+str(js["partition"][s.id])+'}'][1]+=1
                        SeqIO.write(s, ost[l], "fasta")
                    i=i+1
        elif(file_format=="edit"):
            with gzip.open(in_fasta,'rt') as in_handle:
                for line in in_handle:
                    name = line.split()[0]
                    if(name=="root"):
                        for st in ost: #ROOTING=="Origin"
                            st.write(line)
                    elif(js["partition"][name]=="paraphyletic"):
                        para.write(line)
                    else:
                        l = NUMdict[str(js["partition"][name])]
                        if( fasta_count == 0 ):
                            DIRdict['{'+str(js["partition"][name])+'}'][1]+=1
                        ost[l].write(line)
        for st in ost:
            st.close()
        para.close()
    with open(info, 'w') as out:
        out.write(json.dumps(DIRdict))
    for leaf in tree.get_terminals():
        newname=DIRdict[leaf.name][0]
        leaf.name=newname
    Phylo.write(tree, treefile, 'newick')
    return DIRdict

if __name__ == "__main__":
    argvs = sys.argv
    partition_fasta(
        in_fasta_list,
        num_file,
        OUT_DIR,
        wd,
        jpart,
        info,
        treefile,
        subsamplefa,
        ROOTING,
        file_format="fasta"
    )