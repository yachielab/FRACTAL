# TreeAssembly.py
# coded by Naoki Konno

import json
from Bio import Phylo
import sys

color_list=['#ff0000','#ffff00','#00ff00','#00ffff','#0000ff','#ff00ff']

def TreeAssembly(StartDIR, outfname, delete_name):
    init_clade = Phylo.BaseTree.Clade(name=StartDIR)
    tree = Phylo.BaseTree.Tree(init_clade)
    NONTERMINALS=[tree.clade]
    i=0
    while(NONTERMINALS!=[]):
        i+=1
        cstate=NONTERMINALS.pop(0)
        WD=cstate.name
        try:
            downtree=Phylo.read(WD+"/UPSTREAM.nwk",'newick')
            cstate.clades.extend(downtree.clade.clades)
            NONTERMINALS.extend(list(terminal for terminal in downtree.get_terminals()))
        except:
            try:
                downtree=Phylo.read(WD+"/TERMINAL.nwk",'newick')
                if(downtree.clade.clades!=[]):
                    cstate.clades.extend(downtree.clade.clades)
                else:
                    cstate.name=downtree.clade.name
            except:
                print("missing "+WD)
    if (delete_name=="TRUE"):
        for internal_node in tree.get_nonterminals():
            internal_node.name=""
    Phylo.write(tree, outfname, 'newick')

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