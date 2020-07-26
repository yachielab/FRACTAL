from Bio import Phylo
from Bio import SeqIO
import sys
import subprocess
import re
import random
import gzip


def tree_extraction(treefile,nameset, newtreefile):
    tree=Phylo.read(treefile,'newick')
    TipList=tree.get_terminals()
    name_to_numtips={}
    while(TipList!=[]):
        tip=TipList.pop()
        if(tip.name in nameset):
            name_to_numtips[tip.name]=True
        else:
            name_to_numtips[tip.name]=False
    InternalList=tree.get_nonterminals()
    k=0
    while(InternalList!=[]):
        clade=InternalList.pop()
        clade.name="clade"+str(k)
        k+=1
        if(one_of_child_true(clade,name_to_numtips)):
            name_to_numtips[clade.name]=True
        else:
            name_to_numtips[clade.name]=False
    newtree=create_tree(tree, name_to_numtips)
    Phylo.write(newtree, newtreefile, 'newick')
    return newtree

def one_of_child_true(clade,name_to_numtips):
    for child in clade.clades:
        if(name_to_numtips[child.name]):
            return True
    return False

def all_child_true(clade, name_to_numtips):
    for child in clade.clades:
        if(not name_to_numtips[child.name]):
            return False
    return True

def create_tree(tree, name_to_numtips):
    init_clade = Phylo.BaseTree.Clade(name="origin")
    newtree = Phylo.BaseTree.Tree(init_clade)
    stack=[]
    stack.append([tree.clade,init_clade])
    while(stack!=[]):
        statepair=stack.pop()
        cstate=statepair[0]
        parent=statepair[1]
        if(len(cstate.clades)!=1 and all_child_true(cstate,name_to_numtips)):
            parent.clades.append(cstate)
            for child in cstate.clades:
                if(name_to_numtips[child.name]):
                    stack.append([child,cstate])
            cstate.clades=[]
        else:
            for child in cstate.clades:
                if(name_to_numtips[child.name]):
                    stack.append([child,parent])
    return newtree

def fasta_extraction(fastafile,nameset,ext_fastafile):
    with gzip.open(fastafile,'rt') as ihandle, gzip.open(ext_fastafile,'wt') as ohandle:
        sequences=SeqIO.parse(ihandle,'fasta')
        for record in sequences:
            if(record.name in nameset):
                SeqIO.write(record, ohandle,'fasta') 