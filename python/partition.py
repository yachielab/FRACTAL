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
import partition_sequences
import shutil
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

def make_unrooted_after_rooting(treefile,new_treefile,root):
    ref_tree=Phylo.read(treefile,'newick')
    ref_tree.root_with_outgroup(root)
    if len(ref_tree.clade.clades) == 2:
        if len(ref_tree.clade.clades[0].clades) == 2 :
            newtree=Phylo.BaseTree.Tree(ref_tree.clade.clades[0])
            newtree.clade.clades.append(ref_tree.clade.clades[1])
        elif len(ref_tree.clade.clades[1].clades) == 2 :
            newtree=Phylo.BaseTree.Tree(ref_tree.clade.clades[1])
            newtree.clade.clades.append(ref_tree.clade.clades[0])
    Phylo.write(newtree,new_treefile,"newick")

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
    
    # edge2seqlist
    edge2seqlist = {}
    with open(edge_to_sequence_file,'r') as handle:
        for line in handle:
            edge = line.split("\t")[0]
            seq  = line.split("\t")[1].split("\n")[0]
            if (edge in edge2seqlist.keys()):
                edge2seqlist[edge].append(seq)
            else:
                edge2seqlist[edge] = [seq]
    for node in ref_tree.get_terminals() + ref_tree.get_nonterminals():
        if node.name not in edge2seqlist.keys():
            edge2seqlist[node.name] = []
        

    # get paraphyletic group
    for seqname in edge2seqlist[root]:
        partition[seqname]="paraphyletic"
        paraphyletic.append(seqname)
    # get monophyletic groups
    stack.extend(ref_tree.clade.clades[0].clades) #ref_tree.clade.clades[1] : outgroup
    while(stack!=[]):
        cstate=stack.pop()
        d=len(ref_tree.get_path(cstate)) # current depth
        cedge =int(cstate.name.strip('{').strip('}'))
        if(edge2seqlist[cstate.name]==[] and d<depth):
            if(len(cstate.clades)!=0):
                stack.extend(cstate.clades)
            else:
                leaves[cstate.name]=[cstate.name]
        else: # if some seqnames were placed to the edge
            leaves[cstate.name]=list(ter.name for ter in cstate.get_terminals())
            downstream=[]
            downstream.extend(cstate.get_terminals())
            downstream.extend(cstate.get_nonterminals())
            downstream_edge=list(state.name for state in downstream)
            for edge in downstream_edge:
                for seqname in edge2seqlist[edge]:
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

def add_paraphyletic_edit(jpartfname, outputfname, splitted_file_dir, subsample_size, num_of_para):
    # open .jpart file
    with open(jpartfname,"r") as jf:
        jp = jf.read()
    # parse json format
    js = json.loads(jp)

    all_file_pathlist = [ splitted_file_dir + "/" + all_file for all_file in os.listdir(splitted_file_dir) ]
    
    is_gzipped = (all_file_pathlist[0].split(".")[-1] == "gz")
    if (is_gzipped):
        out     = gzip.open(outputfname, 'at')
    else:
        out     = open(outputfname, 'a')
    Nwritten = 0

    for all_file in all_file_pathlist:
        is_gzipped = (all_file.split(".")[-1] == "gz")
        if (is_gzipped):
            ist   = gzip.open(all_file, 'rt')
        else:
            ist   = open(all_file, 'r')
        for line in ist:
            name = line.split()[0]
            if (name != "root"):
                if (len(line.split())>1):
                    editlist_str = line.split()[1]
                else:
                    editlist_str = ""
                if(js["partition"][name]=="paraphyletic"):
                    out.write(line)
                    Nwritten += 1
                    if (Nwritten == num_of_para): break
        ist.close()
        if (Nwritten == num_of_para): break
    out.close()