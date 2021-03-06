'''
Until no paraphyletics, repeat:
    i) repeat 1-4 several times
        1. random sampling (in : FASTA, out : FASTA)
        2. RAxML phylogeny inference (in : FASTA, out : tree)
        3. RAxML parameter optimization (in : FASTA, tree, out RAxML log file)
        4. phylogenetic placement (in : FASTA, Tree, RAxML log file, Query FASTA, out : .jplace file)
    ii) partition.py (in : .jplace files, out : divided FASTAs, directory info file)
'''
import sys
import os
import partition
import subprocess
import shutil
import rename_sequence
import math
import time

def FRACluster(WD, MAX_ITERATION, SUBSAMPLE_SIZE, NODESDIR, THRESHOLD, THREAD_NUM, NUMFILE, QSUBDIR, CODEDIR, ROOTING, MODEL, OPTION,TREEMETHOD, ALIGNMETHOD, EPANG, RAXMLSEQ, RAXMLPAR, SOFTWARE,NODE_COUNT,INIT_SEQ_COUNT,SEED):
    start=time.time() # in order to get the time which one cycle takes
    subprocess.call("which bash",shell=True)
    os.chdir(WD) # move to Working Directory
    # check if outgroup exists or not (sequence named "root")
    if rename_sequence.outgroup_check_fast("INPUT.fa"):
        print("INPUT.fa properly include outgroup sequence.")
    else:
        print("No sequence named \"root\" in INPUT.fa! INPUT.fa should include outgroup sequence named \"root\".")
        return
    
    # check number of sequences
    seq_array=rename_sequence.count_sequence_fast("INPUT.fa")
    seq_count=seq_array[0]
    if(INIT_SEQ_COUNT==0): INIT_SEQ_COUNT=seq_count # only in d0
    seq_length=seq_array[1]
    raxml_thread_num=min(max(seq_length//500,2),THREAD_NUM) # use 1 thread per 500 bp in RAxML
    depth=max(math.floor(math.log2(seq_count/THRESHOLD))+2,2)
    tree_thread_num=THREAD_NUM
    if(TREEMETHOD=="raxmlML" or TREEMETHOD=="raxmlMP"): tree_thread_num=raxml_thread_num

    i=1
    # call direct tree reconstruction
    if(seq_count<=THRESHOLD):
        if(seq_count<4):
            partition.tiny_tree("INPUT.fa","TERMINAL.nwk")
            print("seq_count < 4!")
        else:
            os.mkdir("TREE")
            os.mkdir("PARAM")
            os.chdir(WD+"/TREE")
            subprocess.call("bash "+CODEDIR+"/shell/TREE.sh -n "+str(tree_thread_num)+" -m "+TREEMETHOD+" -a "+ALIGNMETHOD+" -f "+WD+"/INPUT.fa -c "+CODEDIR+" -w "+WD+"/TREE -p \""+str(OPTION)+"\" -d "+MODEL+" -q "+SOFTWARE,shell=True)
            partition.rooting_and_remove(WD+"/INPUT.fa.aligned.tree",WD+"/TERMINAL.nwk","root")

    # call fractal FRACTAL
    elif(NODE_COUNT>1 and seq_count<=INIT_SEQ_COUNT//NODE_COUNT): # quit distribution after the available computer node saturated 
        os.mkdir("TREE")
        os.mkdir("PARAM")
        os.chdir(WD+"/TREE")
        if (TREEMETHOD!="unspecified"):
            subprocess.call("bash "+CODEDIR+"/FRACTAL.sh "+"-x "+str(MAX_ITERATION)+" -k "+str(SUBSAMPLE_SIZE)+" -i "+WD+"/INPUT.fa" +" -t "+str(THRESHOLD)+" -m "+TREEMETHOD+" -b "+MODEL+" -c "+str(THREAD_NUM)+" -e",shell=True)
        else:
            subprocess.call("bash "+CODEDIR+"/FRACTAL.sh "+"-x "+str(MAX_ITERATION)+" -k "+str(SUBSAMPLE_SIZE)+" -i "+WD+"/INPUT.fa" +" -t "+str(THRESHOLD)+" -s "+SOFTWARE+" -b "+MODEL+" -c "+str(THREAD_NUM)+" -e",shell=True)
        shutil.move(WD+"/TREE/FRACTALout.nwk",WD+"/TERMINAL.nwk")

    # call FRACTAL cycle
    else:
        os.mkdir("SUBSAMPLE")
        os.mkdir("TREE")
        os.mkdir("PARAM")
        os.mkdir("EPANG")
        os.mkdir("PARTITION")
        os.mkdir("ANCSEQ")
        prev_para=seq_count

        while i<MAX_ITERATION:
            #################
            #random sampling#
            #################
            if(os.path.isfile("ITERATION.fa")):
                rename_sequence.random_sampling("ITERATION.fa","SUBSAMPLE/SUBSAMPLE.fa",SUBSAMPLE_SIZE,seed=SEED)
                shutil.rmtree("EPANG")
                os.mkdir("EPANG")
                shutil.rmtree("TREE")
                os.mkdir("TREE")
            else:
                rename_sequence.random_sampling("INPUT.fa","SUBSAMPLE/SUBSAMPLE.fa",SUBSAMPLE_SIZE,seed=SEED,n = seq_count)
            #################
            #rename sequence#
            #################
            rename_sequence.rename_sequence("SUBSAMPLE/SUBSAMPLE.fa", "SUBSAMPLE/RENAMED_"+str(i)+".fa")
            #######################################
            #construct subsample tree as reference#
            #######################################
            os.chdir(WD+"/TREE")
            subprocess.call("bash "+CODEDIR+"/shell/TREE.sh -n "+str(tree_thread_num)+" -m "+TREEMETHOD+" -a "+ALIGNMETHOD+" -f "+WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa -c "+CODEDIR+" -w "+WD+"/TREE -p \""+str(OPTION)+"\" -d "+MODEL+" -q "+SOFTWARE,shell=True)
            os.chdir(WD+"/PARAM")
            if(seq_count>=30):
                subprocess.call(RAXMLPAR+" -T "+str(raxml_thread_num)+" -f e -s "+WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned -t "+WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned.tree -n PARAM_"+str(i)+" -m " +MODEL,shell=True)
            else:
                subprocess.call(RAXMLSEQ+" -f e -s "+WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned -t "+WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned.tree -n PARAM_"+str(i)+" -m " +MODEL,shell=True)
            ########################################
            #Phylogenetic placement & visualization#
            ########################################
            os.chdir(WD)
            if(seq_count>10000): nodenum = (NODE_COUNT*seq_count)//INIT_SEQ_COUNT-1
            else: nodenum = 0
            rename_sequence.distributed_EPAng(WD, EPANG, WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned", WD+"/PARAM/RAxML_result.PARAM_"+str(i), WD+"/PARAM/RAxML_info.PARAM_"+str(i), WD+"/INPUT.fa", WD+"/EPANG", THREAD_NUM, nodenum,CODEDIR,seq_count)
            ####################
            #parse .jplace file#
            ####################
            os.chdir(WD)
            para=partition.partition(WD+"/EPANG/placement_tree.out",WD+"/EPANG/edge_to_seqname_all.out",WD+"/PARTITION/partition"+str(i)+".out",depth)
            print("detected "+str(para)+" paraphyletic sequences")
            ##################################################
            #get paraphyletic sequences and make ITERATION.fa#
            ##################################################
            if(para>prev_para):
                i-=1 # previous subsample achieved minimum number of paraphyletic groups
                break
            if(para!=0):
                shutil.copyfile("SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned","ITERATION.fa")
                partition.add_paraphyletic_fa(WD+"/PARTITION/partition"+str(i)+".out","ITERATION.fa", "INPUT.fa",SUBSAMPLE_SIZE,para)
                i+=1
                prev_para=para
            else:
                break
        ##################
        #partition .fasta#
        ##################
        os.chdir(WD)
        DIRdict=partition.partition_fasta("INPUT.fa",NUMFILE,NODESDIR,WD,WD+"/PARTITION/partition"+str(min(i,MAX_ITERATION-1))+".out","PARTITION.info","UPSTREAM.nwk",WD+"/ANCSEQ/RAxML_marginalAncestralStates.ANCSEQ", WD+"/ANCSEQ/RAxML_nodeLabelledRootedTree.ANCSEQ", WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa",ROOTING)
        partition.qsub_prep(MAX_ITERATION, SUBSAMPLE_SIZE, NODESDIR, THRESHOLD, THREAD_NUM, NUMFILE, QSUBDIR, CODEDIR, DIRdict,ROOTING,MODEL,OPTION,TREEMETHOD, ALIGNMETHOD,EPANG, RAXMLSEQ, RAXMLPAR,SOFTWARE,NODE_COUNT,INIT_SEQ_COUNT,SEED)
        ##################
        #delete files    #
        ##################
        os.remove("INPUT.fa")
        os.remove("tmp.nwk")
        shutil.rmtree("ANCSEQ")
        shutil.rmtree("EPANG")
        shutil.rmtree("PARAM")
        shutil.rmtree("PARTITION")
        shutil.rmtree("SUBSAMPLE")
        shutil.rmtree("TREE")
    
    elapsed_time=time.time()-start
    with open(WD+"/time.out", 'w') as handle:
        handle.write(str(seq_count)+","+str(elapsed_time) + ",sec,"+str(i)+",subsamplings\n")

if __name__ == "__main__":
    argvs = sys.argv
    FRACluster(argvs[1],int(argvs[2]),int(argvs[3]),argvs[4],int(argvs[5]),int(argvs[6]),argvs[7],argvs[8], argvs[9],argvs[10],argvs[11],argvs[12],argvs[13],argvs[14],argvs[15],argvs[16],argvs[17],argvs[18],int(argvs[19]),int(argvs[20]),argvs[21])