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

def add_paraphyletic_fa(jpartfname, outputfname, aligned_all_dir, subsample_size, num_of_para, file_format = "fasta"):
    # open .jpart file
    with open(jpartfname,"r") as jf:
        jp = jf.read()
    # parse json format
    js = json.loads(jp)

    all_fa_pathlist = [ aligned_all_dir + "/" + all_fa for all_fa in os.listdir(aligned_all_dir) ]
    for all_fa in all_fa_pathlist:
        is_gzipped = (all_fa.split(".")[-1] == "gz")
        if (is_gzipped):
            allfa   = gzip.open(all_fa, 'rt')
            out     = gzip.open(outputfname, 'at')
        else:
            allfa   = open(all_fa, 'r')
            out     = open(outputfname, 'a')
        # add paraphyletic sequences into subsample
        if (file_format == "fasta"):
            Nwritten = 0
            handle = SeqIO.parse(allfa, "fasta")
            for record in handle:                                    # should be faster!
                if (Nwritten == num_of_para):
                    break
                if(record.id!="root"):
                    if(js["partition"][record.id]=="paraphyletic"): 
                        SeqIO.write(record, out, "fasta")
                        Nwritten += 1
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
        allfa.close()
        out.close()

def get_ancseq(ancseq,ancnum):
    ancname=str(ancnum)
    with open(ancseq, 'r') as handle:
        for line in handle:
            ls=line.split()
            if(ls[0]==ancname):
                return ls[1]
    print("no sequence named "+ancname)

def partition_fasta(
    in_fasta_dirlist,
    num_file,
    OUT_DIR,
    wd,
    jpart,
    treefile,
    ROOTING,
    file_format="fasta", 
    nodenum=1, 
    codedir=None
    ):
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
    
    dirpath_set = set()
    with open(wd + "/seqname_dirpath.txt", 'w') as handle:
        for seqname in list(js["partition"].keys()):
            if js["partition"][seqname] == "paraphyletic":
                handle.write(seqname +'\t' + wd + '\tparaphyletic\n')
            else:
                dirpath = DIRdict['{'+str(js["partition"][seqname])+'}'][0]
                dirpath_set.add(dirpath)
                handle.write(seqname +'\t' + dirpath + '\n')
    subprocess.call(
        "cat " + wd + "/seqname_dirpath.txt | cut -f2 | sort | uniq -c > " + wd + "/Nseq_dirpath.txt",
        shell = True
    )

    for fasta_count, splitted_fasta_dir in enumerate(in_fasta_dirlist):

        is_gzipped = (os.listdir(splitted_fasta_dir)[0].split(".")[-1] == "gz")
        splitted_fasta_list = os.listdir(splitted_fasta_dir)
        splitted_fpath_list = [splitted_fasta_dir + "/" + splitted_fasta for splitted_fasta in splitted_fasta_list]
        Nfiles_total = len(splitted_fasta_list)
        dirpath_list = list(sorted([wd] + list(dirpath_set)))

        if (is_gzipped):
            gzip_extention = ".gz"
            gzip_command   = "|gzip"
        else:
            gzip_extention = ""
            gzip_command   = ""

        if (file_format=="fasta"):

            if (nodenum > 1):
                
                # assign splitted files to each node: same as distributed placement
                #Nfiles_per_node = len(splitted_fasta_list) // nodenum # Only the last node may treat more number of files
                node2filelist = []
                for i in range(nodenum):
                    node2filelist.append([])
                for j, file_name in enumerate(splitted_fasta_list):
                    node2filelist[j%nodenum].append(file_name)
                print ("node2filelist", node2filelist)

                dname = wd.split("/")[-1]
                for i in range(nodenum):
                    with open(wd+"/../../qsub_dir/qsub_"+dname+"."+str(i)+".partition.sh", 'w') as handle:
                        PATH = (subprocess.\
                                    Popen(
                                        'echo $PATH',
                                        stdout=subprocess.PIPE,
                                        shell=True
                                    ).communicate()[0]
                                ).decode('utf-8')
                        PATH = (PATH.split('\n'))[0]
                        LD_LIBRARY_PATH = (
                            subprocess.\
                                Popen(  
                                    'echo $LD_LIBRARY_PATH', 
                                    stdout=subprocess.PIPE,
                                    shell=True
                                ).communicate()[0]
                            ).decode('utf-8')
                        LD_LIBRARY_PATH = (LD_LIBRARY_PATH.split('\n'))[0]
                        handle.write("#!/bin/bash\n")
                        handle.write("#$ -S /bin/bash\n")
                        handle.write("PATH={}\n".format(PATH))
                        handle.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
                        handle.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
                        inputFASTA_filepathlist = [splitted_fasta_dir + "/" + splitted_file for splitted_file in node2filelist[i] ]
                        handle.write(
                            "python3 "                               +
                            codedir + "/python/partition_sequences.py "  +
                            ":".join(inputFASTA_filepathlist) + " "  +
                            ":".join(dirpath_list)+" "               +
                            wd + "/seqname_dirpath.txt"              +
                            "\n"
                            )
                        handle.write(
                            "echo finished > " + splitted_fasta_dir+"/"+dname+"."+str(i)+".partition.sh.finished\n"
                            )
                # wait for all partition jobs finish
                Nunclassified    = len(os.listdir(splitted_fasta_dir))
                Nfinished        = 0
                while(Nunclassified != 0):
                    file_list= os.listdir(splitted_fasta_dir)
                    Nunclassified    = 0
                    Nfinished        = 0
                    for filename in file_list:
                        if (filename.split(".")[-1] == 'finished'):
                            Nfinished += 1
                        elif (filename.split(".")[-1] == 'fa' or filename.split(".")[-1]=='gz'):
                            Nunclassified += 1
                    if (Nunclassified > 0 and Nfinished == nodenum):
                        subprocess.call(
                            "mv "+wd+"/../../executed/qsub_"+dname+".*.partition.sh " + wd+"/../../qsub_dir",
                            shell=True
                        )
                        
            else: # sequential mode
                print(splitted_fpath_list, dirpath_list, wd + "/seqname_dirpath.txt")
                partition_sequences.partition_sequences(splitted_fpath_list, dirpath_list, wd + "/seqname_dirpath.txt")
            for dirpath in dirpath_list:
                if (dirpath != wd):
                    subprocess.call(
                        "cat " + dirpath + "/*.count > "+dirpath+"/file2Nseq.txt; "+
                        "rm "+dirpath + "/*.count; "                               +
                        "cat "+wd+"/root.fa > "+dirpath+"/root.fa",
                        shell=True
                    )
            problematic_filenames      = wd + "/*.part*fa"+gzip_extention
            problematic_concatfilename = wd+"/INPUT.fa.problematic"+gzip_extention
            subprocess.call(
                "cat "                                   +
                problematic_filenames                    +
                "|seqkit grep -r -p ^root -v"            +
                ">> " + problematic_concatfilename + ";" +
                "rm " + problematic_filenames            +
                "",
                shell = True
                )

        elif(file_format=="edit"): # TO DO
            ost=[]
            for i in range(num_mono):
                ost.append(gzip.open(OUT_DIR+"/d"+str(num+i)+"/"+in_fasta.split("/")[-1],'wt'))
            para=gzip.open(wd+"/"+in_fasta.split("/")[-1]+".problematic.gz",'wt')

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
    #with open(info, 'w') as out:
    #    out.write(json.dumps(DIRdict))
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
                ARGVS[28] = seq_count_when_aligned


            command="python3 "
            for arg in ARGVS: 
                if str(arg)=="":
                    arg="\"\""
                command += str(arg) + " "
            qf.write(command)

def tiny_tree(INPUTfile,OUTPUTnwk, file_format="fasta"):
    is_gzipped = (INPUTfile.split(".")[-1] == "gz")
    if is_gzipped:
        handle  = gzip.open(INPUTfile, 'rt')
    else:
        handle  = open(INPUTfile, 'r')

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

    handle.close()