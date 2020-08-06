# jplace_parse.py
# coded by Naoki Konno
# parse .jplace file and returns dictionary of 'tree', 'placement dictionary', and 'placement list'

import json
from Bio import Phylo
from io import StringIO
import sys
import re
import random

def correspond(treestr):
    corr={}
    pattern1 = r"s\d+"
    pattern2 = r"{\d+}"
    i=0
    root=''
    while(i<len(treestr)):
        matchOB=re.search(pattern1, treestr[i:])
        if matchOB:
            name=matchOB.group()
            i+=matchOB.end()+1
            matchOB=re.search(pattern2, treestr[i:])
            if matchOB:
                num=matchOB.group()
                corr[num]=name
                if name=='s0':
                    root=num
                i+=matchOB.end()+1
            else:
                print("ERR: no corresponding number of "+name)
                break
        else:
            break
    return [corr,root]

def parse_jplace(fname, placement_method,seed):
    if(len(seed)!=0):random.seed(int(seed))
    with open(fname,"r") as jf:
        jp = jf.read()
    # parse json format
    js = json.loads(jp)
    # get tree information
    treestr = js['tree']
    tree=Phylo.read(StringIO(treestr),'newick')
    # get placement information
    jdict = js['placements']
    placement_list = [] # edge name -> sequence name
    for i in range(len(tree.get_terminals())+len(tree.get_nonterminals())):
        placement_list.append([])
    for placement in jdict:
        if(placement_method=="epa-ng"):
            problist = list(pl[2] for pl in placement['p'])
            maxidx=problist.index(max(problist))
            edge = placement['p'][maxidx][0]
            name = placement['n'][0]
        elif(placement_method=="epa_MP"):
            tree.clade.name = "tree_top"
            equally_parsimonious_edge_list = list('{'+str(pl[0])+'}' for pl in placement['p'])
            #edge = random.choice(equally_parsimonious_edge_list)
            edge_str = tree.common_ancestor(equally_parsimonious_edge_list).name
            print(equally_parsimonious_edge_list, edge_str)
            if (edge_str == "tree_top"):
                root=correspond(treestr)[1]
                edge_str = '{' + str(root) + '}'
            edge = int(edge_str.split('{')[1].split('}')[0])
            name = placement['n'][0]
        if(name!='root'):
            placement_list[edge].append(name)
    with open('placement_tree.out','w') as handle:
        handle.write(treestr)
    with open("edge_to_seqname.out",'w') as handle:
        for seqnamelist in placement_list:
            for seqname in seqnamelist:
                handle.write(seqname+",")
            handle.write("\n")

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    #main function
    argvs=sys.argv
    parse_jplace(argvs[1],argvs[2],argvs[3])