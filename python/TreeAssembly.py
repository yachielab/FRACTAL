# TreeAssembly.py
# coded by Naoki Konno

import json
from Bio import Phylo
import sys
import os
import extraction

color_list=['#ff0000','#ffff00','#00ff00','#00ffff','#0000ff','#ff00ff']

def TreeAssembly(StartDIR, outfname, delete_name):
    init_clade = Phylo.BaseTree.Clade(name=StartDIR)
    tree = Phylo.BaseTree.Tree(init_clade)
    NONTERMINALS=[tree.clade]
    i=0

    print("Start tree assembly...")

    tip_name_list = []
    extraction_is_needed = False
    while(NONTERMINALS!=[]):
        i+=1
        cstate=NONTERMINALS.pop(0)
        WD=cstate.name
        if (os.path.exists(WD+"/UPSTREAM.nwk")):
            downtree=Phylo.read(WD+"/UPSTREAM.nwk",'newick')
            cstate.clades.extend(downtree.clade.clades)
            NONTERMINALS.extend(list(terminal for terminal in downtree.get_terminals()))
        elif (os.path.exists(WD+"/TERMINAL.nwk")):
            downtree=Phylo.read(WD+"/TERMINAL.nwk",'newick')
            tip_name_list.extend([tip.name for tip in downtree.get_terminals()])

            if(downtree.clade.clades!=[]):
                cstate.clades.extend(downtree.clade.clades)
            else:
                cstate.name=downtree.clade.name
        else:
            print("No sequences in clade", cstate.name)
            extraction_is_needed = True

    if (delete_name=="TRUE"):
        for internal_node in tree.get_nonterminals():
            internal_node.name=""

    print("Finished tree assembly")

    if (extraction_is_needed):
        print("Start tree extraction...")
        
        tree = extraction.tree_extraction_biopython(tree,set(tip_name_list))

        print("Finished tree extraction")
    
    print("Start writing a newick file...")

    Phylo.write(tree, outfname, 'newick')

    print("Finish writing a newick file")

'''
main function
        commandline argument: "python3 TreeAssembly.py <StartDIR> <outfname> <coloringfilepath>"
'''

if __name__ == "__main__":

    try:
        sys.setrecursionlimit(1000000000)
    except:
        None

    argvs = sys.argv
    TreeAssembly(argvs[1],argvs[2],argvs[3])