# jplace_parse.py
# coded by Naoki Konno
# parse .jplace file and returns dictionary of 'tree', 'placement dictionary', and 'placement list'

import json
from Bio import Phylo, SeqIO
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

def parse_jplace(fname, placement_method, infasta_fpath, seed, careful=1):
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

    problematic_set = set()
    for placement in jdict:
        if(placement_method=="epa-ng"):
            edge_prob_list = [{'edge':pl[0], 'prob':pl[2]} for pl in placement['p']]
            if (careful>1 and len(edge_prob_list)>1 ):
                tree.clade.name = "tree_top"
                root=correspond(treestr)[1]
                edge_prob_list_sorted = sorted(edge_prob_list, key=lambda x:x['prob'], reverse=True)
                best_places  = ['{'+str(edge_prob['edge'])+'}' for edge_prob in edge_prob_list_sorted[0:careful]]
                edge_str = tree.common_ancestor(best_places).name
                if (edge_str == "tree_top"):
                    edge_str = root
                edge = int(edge_str.split('{')[1].split('}')[0])
            else:
                problist = list(pl[2] for pl in placement['p'])
                maxidx=problist.index(max(problist))
                edge = placement['p'][maxidx][0]
            name = placement['n'][0]
        elif(placement_method=="epa_MP"):
            tree.clade.name = "tree_top"
            equally_parsimonious_edge_list = list('{'+str(pl[0])+'}' for pl in placement['p'])
            if (careful==1):
                edge_str = random.choice(equally_parsimonious_edge_list)
            else:
                edge_str = tree.common_ancestor(equally_parsimonious_edge_list).name
                if (edge_str == "tree_top"):
                    root=correspond(treestr)[1]
                    edge_str = root
            edge = int(edge_str.split('{')[1].split('}')[0])
            name = placement['n'][0]
        if(name!='root'):
            placement_list[edge].append(name)

        if('{'+str(edge)+'}' == root):
            problematic_set.add(name)

    with open('placement_tree.out','w') as handle:
        handle.write(treestr)
    '''
    with open("edge_to_seqname.out",'w') as handle:
        for seqnamelist in placement_list:
            for seqname in seqnamelist:
                handle.write(seqname+",")
            handle.write("\n")
    '''
    with open("edge_to_seqname.out",'w') as handle:
        for i, seqnamelist in enumerate(placement_list):
            for seqname in seqnamelist:
                handle.write('{'+str(i)+'}'+"\t"+seqname+"\n")

    with open("problematic.fa", 'w') as handle:
        records = SeqIO.parse(infasta_fpath,'fasta')
        for record in records:
            if(record.name in problematic_set):
                SeqIO.write(record, handle, 'fasta')

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    #main function
    argvs=sys.argv
    parse_jplace(argvs[1],argvs[2],argvs[3],argvs[4])