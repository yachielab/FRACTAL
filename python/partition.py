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
#from memory_profiler import profile

def rooting(nwkfilepath,newnwkpath,root):
    with open(nwkfilepath) as handle:
        tree1=Phylo.read(handle,"newick")
        tree1.root_with_outgroup(root)
        tree2=Phylo.BaseTree.Tree(tree1.clade)
    Phylo.write(tree2,newnwkpath,"newick")

def rooting_and_remove(nwkfilepath,newnwkpath,root):
    with open(nwkfilepath) as handle:
        tree1=Phylo.read(handle,"newick")
        tree1.root_with_outgroup(root)
        tree2=Phylo.BaseTree.Tree(tree1.clade[0])
    Phylo.write(tree2,newnwkpath,"newick")

def make_unrooted(treefile,new_treefile):
    ref_tree=Phylo.read(treefile,'newick')
    if len(ref_tree.clade.clades) == 2:
        if len(ref_tree.clade.clades[0].clades) == 2 :
            ref_tree.clade = ref_tree.clade.clades[0]
            ref_tree.clade.clades.append(ref_tree.clade.clades[1])
        elif len(ref_tree.clade.clades[1].clades) == 2 :
            ref_tree.clade = ref_tree.clade.clades[1]
            ref_tree.clade.clades.append(ref_tree.clade.clades[0])
    Phylo.write(ref_tree,new_treefile,"newick")

def partition(treefile, edge_to_sequence_file, jpartitionfname, depth):
    partition={}; leaf_to_Nseq={}
    paraphyletic=[]
    stack=[]
    leaves={}
    
    # get tree information
    with open(treefile,'r') as handle:
        treestr = handle.read()
    # get relation of 'sXXX' <--> {YYY}
    relation=jplace_parse.correspond(treestr)
    corr=relation[0]
    root=relation[1]
    ref_tree=Phylo.read(treefile,'newick')
    ref_tree.root_with_outgroup(root)
    print(ref_tree)
    
    # place list
    place_list=[]
    with open(edge_to_sequence_file,'r') as handle:
        for line in handle:
            line=line.split("\n")[0]
            if(len(line)>0): line=line[0:(len(line)-1)];place_list.append(line.split(","))
            else: place_list.append([])

    # get paraphyletic group
    for seqname in place_list[int(root.strip('{').strip('}'))]:
        partition[seqname]="paraphyletic"
        paraphyletic.append(seqname)
    # get monophyletic groups
    stack.extend(ref_tree.clade.clades[0].clades) #ref_tree.clade.clades[1] : outgroup
    while(stack!=[]):
        cstate=stack.pop()
        d=len(ref_tree.get_path(cstate)) # current depth
        if (cstate.name == None):
            if(d<depth):
                if(len(cstate.clades)!=0):
                    stack.extend(cstate.clades)
                else:
                    leaves[cstate.name]=[cstate.name]
        else:
            cedge =int(cstate.name.strip('{').strip('}'))
            if(place_list[cedge]==[] and d<depth):
                if(len(cstate.clades)!=0):
                    stack.extend(cstate.clades)
                else:
                    leaves[cstate.name]=[cstate.name]
            else: # if some seqnames were placed to the edge
                leaves[cstate.name]=list(ter.name for ter in cstate.get_terminals())
                downstream=[]
                downstream.extend(cstate.get_terminals())
                downstream.extend(cstate.get_nonterminals())
                downstream_edge=list(int(state.name.strip('{').strip('}')) for state in downstream)
                for edge in downstream_edge:
                    for seqname in place_list[edge]:
                        if seqname != "root":
                            partition[seqname]=cedge
                            if(cedge in leaf_to_Nseq.keys()): leaf_to_Nseq[cedge]+=1
                            else: leaf_to_Nseq[cedge]=1
                cstate.clades=[]
    tree=Phylo.BaseTree.Tree(ref_tree.clade.clades[0])
    Phylo.write(tree, 'tmp.nwk', 'newick')
    with open('tmp.nwk','r') as tmp:
        treestr=tmp.read()
    outputdict = {}
    outputdict["tree"]=treestr
    outputdict["paraphyletic"]=paraphyletic
    outputdict["partition"]=partition
    outputdict["leaves"]=leaves
    outputdict["corr"]=corr
    with open(jpartitionfname,'w') as out:
        out.write(json.dumps(outputdict))
    return len(paraphyletic), max(list(leaf_to_Nseq.values()))

def add_paraphyletic_fa(jpartfname, outputfname, all_fa, subsample_size, num_of_para, file_format = "fasta"):
    # open .jpart file
    with open(jpartfname,"r") as jf:
        jp = jf.read()
    # parse json format
    js = json.loads(jp)
    '''
    # generate random indices
    if(subsample_size<num_of_para):
        rand_idx=random.sample(range(num_of_para),subsample_size)
    else:
        rand_idx=range(num_of_para)
    '''
    # add paraphyletic sequences into subsample
    with open(all_fa,'r') as allfa, open(outputfname,'a') as out:
        if (file_format == "fasta"):
            handle = SeqIO.parse(allfa, "fasta")
            for record in handle:
                if(record.id!="root"):
                    if(js["partition"][record.id]=="paraphyletic"):
                        SeqIO.write(record, out, "fasta")
        elif (file_format=="edit"):
            for line in allfa:
                name         = line.split()[0]
                if (name != "root"):
                    if (len(line.split())>1):
                        editlist_str = line.split()[1]
                    else:
                        editlist_str = ""
                    if(js["partition"][name]=="paraphyletic"):
                        out.write(line)

def get_ancseq(ancseq,ancnum):
    ancname=str(ancnum)
    with open(ancseq, 'r') as handle:
        for line in handle:
            ls=line.split()
            if(ls[0]==ancname):
                return ls[1]
    print("no sequence named "+ancname)

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
            ost.append(open(OUT_DIR+"/d"+str(num+i)+"/"+in_fasta.split("/")[-1],'w'))
        para=open(wd+"/"+in_fasta.split("/")[-1]+".problematic",'w')

        if (file_format=="fasta"):
            with open(in_fasta,'r') as in_handle:
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
            with open(in_fasta,'r') as in_handle:
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

def qsub_prep(ARGVS, QSUBDIR, DIRdict, INITIAL_SEQ_COUNT, seq_count_when_aligned):
    for key in DIRdict.keys():
        ls=DIRdict[key][0].split("/")
        num=ls[len(ls)-1]
        
        PATH = (
            subprocess.Popen(
                'echo $PATH',
                stdout=subprocess.PIPE,
                shell=True)
                .communicate()[0]
                ).decode('utf-8')

        PATH = PATH.split('\n')[0]
        
        LD_LIBRARY_PATH = (
            subprocess.Popen(
                'echo $LD_LIBRARY_PATH',
                stdout=subprocess.PIPE,
                shell=True
                ).communicate()[0]
                ).decode('utf-8')
        
        LD_LIBRARY_PATH = (LD_LIBRARY_PATH.split('\n'))[0]
        
        with open(QSUBDIR+"/qsub_"+num+".sh", 'w') as qf:
            qf.write("#!/bin/bash\n")
            qf.write("#$ -S /bin/bash\n")
            qf.write("PATH={}\n".format(PATH))
            qf.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
            
            # change arguments for the next FRACTAL iteration
            ARGVS[1]  = DIRdict[key][0]
            ARGVS[20] = INITIAL_SEQ_COUNT
            if (ARGVS[14] == "unaligned"):
                ARGVS[26] = seq_count_when_aligned


            command="python3 "
            for arg in ARGVS: 
                if str(arg)=="":
                    arg="\"\""
                command += str(arg) + " "
            qf.write(command)

def tiny_tree(INPUTfile,OUTPUTnwk, file_format="fasta"):
    with open(INPUTfile,'r') as handle:
        names=[]
        if (file_format == "fasta"):
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                if(record.id!="root"):
                    names.append(record.id)
        elif(file_format == "edit"):
            for line in handle:
                name = line.split()[0]
                if (name != "root"): names.append(name)

    if(len(names)==1):
        init_clade = Phylo.BaseTree.Clade(name=names[0])
        tree = Phylo.BaseTree.Tree(init_clade)
    elif(len(names)==2):
        init_clade = Phylo.BaseTree.Clade()
        tree = Phylo.BaseTree.Tree(init_clade)
        tree.clade.clades.extend(list(Phylo.BaseTree.Clade(name=name) for name in names))
    else:
        print("tiny_tree() Error : len(names)=")
        print(len(names))
    Phylo.write(tree, OUTPUTnwk, 'newick')