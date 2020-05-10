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
import placement
import math
import time

def FRACluster(ARGVS, WD, MAX_ITERATION, SUBSAMPLE_SIZE, NODESDIR, THRESHOLD, THREAD_NUM, NUMFILE, QSUBDIR, CODEDIR, 
               ROOTING, MODEL, OPTION,TREEMETHOD, ALIGNED, EPANG, RAXMLSEQ, RAXMLPAR, SOFTWARE,NODE_COUNT,
               INIT_SEQ_COUNT,SEED,ML_or_MP, 
               ALIGNER="unspecified", HMM_PROFILER="unspecified", HMM_ALIGNER="unspecified",
               seq_count_when_aligned=None
               ):
    
    ######## parameter ########
    ALIGNMENT_TIMING_PARAMETER = 0.5
    ###########################

    start = time.time() # in order to get the time which one cycle takes
    
    os.chdir(WD) # move to Working Directory
    
    # check if outgroup exists or not (sequence named "root")
    if rename_sequence.outgroup_check_fast(WD+"/INPUT.fa"):
        print("INPUT.fa properly include outgroup sequence")
    else:
        print("No sequence named \"root\" in INPUT.fa!\
               INPUT.fa should include outgroup sequence named \"root\"")
        return
    
    # check number of sequences
    seq_array                 = rename_sequence.count_sequence_fast("INPUT.fa")
    seq_count                 = seq_array[0]
    if(INIT_SEQ_COUNT==0): 
        INIT_SEQ_COUNT        = seq_count # only in d0
        seq_count_when_aligned= None
    seq_length                = seq_array[1]
    raxml_thread_num          = min(max(seq_length//500,2),THREAD_NUM) # use 1 thread per 500 bp in RAxML
    depth                     = max(math.floor(math.log2(seq_count/THRESHOLD))+2,2)
    tree_thread_num           = THREAD_NUM
    if(TREEMETHOD=="raxmlML" or
       TREEMETHOD=="raxmlMP"    ): 
        tree_thread_num       = raxml_thread_num

    # check if aligned
    if (ALIGNED=="unaligned"):
        if(os.path.isfile(WD+"/INPUT.fa.aligned")):
            if (seq_count < seq_count_when_aligned * ALIGNMENT_TIMING_PARAMETER):

                print(WD,seq_count,seq_count_when_aligned,"alignment needed!!")

                shutil.remove(WD+"/INPUT.fa.aligned")
                seq_count_when_aligned = seq_count
                
            else:
                INPUT_FA = WD+"/INPUT.fa.aligned"
                ALIGNED  = "aligned"
        else:

            print(WD,seq_count,seq_count_when_aligned,"alignment needed!!")

            INPUT_FA               = WD+"/INPUT.fa"
            seq_count_when_aligned = seq_count
    i=1
    
    # call direct tree reconstruction

    if(seq_count<=THRESHOLD):
        if(seq_count<4):
            partition.tiny_tree("INPUT.fa","TERMINAL.nwk")
            print("seq_count < 4!")
        else:
            os.mkdir(    "TREE")
            os.chdir(WD+"/TREE")
            subprocess.call("bash " + CODEDIR + "/shell/TREE.sh" +
                            " -n "  + str(tree_thread_num)       +
                            " -m "  + TREEMETHOD                 +
                            " -a "  + ALIGNED                    +
                            " -f "  + WD+"/INPUT.fa"             +
                            " -c "  + CODEDIR                    +
                            " -w "  + WD+"/TREE"                 +
                            " -p \""+ str(OPTION)+"\""           +
                            " -d "  + MODEL                      +
                            " -q "  + SOFTWARE                   +
                            " -s "  + ALIGNER                    ,
                            shell   = True                       )
            partition.rooting_and_remove(WD+"/INPUT.fa.aligned.tree",WD+"/TERMINAL.nwk","root")

    # call fractal FRACTAL # don't forget to change!!!!!
    elif(NODE_COUNT>1 and seq_count<=INIT_SEQ_COUNT//NODE_COUNT): 
        # Quit distribution after the available computer node saturated 
        os.mkdir("TREE")
        os.chdir(WD+"/TREE")
        FRACTAL_COMMAND = "bash "+ CODEDIR + "/FRACTAL.sh"  + \
                          " -x " + str(MAX_ITERATION)       + \
                          " -k " + str(SUBSAMPLE_SIZE)      + \
                          " -i " + WD + "/INPUT.fa"         + \
                          " -t " + str(THRESHOLD)           + \
                          " -b " + MODEL                    + \
                          " -c " + str(THREAD_NUM)          + \
                          " -p " + ML_or_MP                 + \
                          " -a " + OPTION                   + \
                          " -r " + SEED
        
        if (TREEMETHOD!="unspecified"): 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -m "+TREEMETHOD
        else: 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -s "+SOFTWARE
        if (ALIGNED=='unaligned'):
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -u "
        
        subprocess.call(FRACTAL_COMMAND,shell=True)
        shutil.move(WD+"/TREE/FRACTALout.nwk",WD+"/TERMINAL.nwk")

    # call FRACTAL cycle
    else:
        os.mkdir("SUBSAMPLE")
        os.mkdir("TREE")
        os.mkdir("PARAM")
        os.mkdir("EPANG")
        os.mkdir("PARTITION")
        os.mkdir("ANCSEQ")
        
        prev_para = seq_count

        while i<MAX_ITERATION:
            
            #################
            #random sampling#
            #################
            if(os.path.isfile(WD+"/INPUT.fa.aligned")):
                INPUT_FA = WD+"/INPUT.fa.aligned"

            if(os.path.isfile("ITERATION.fa")):
                rename_sequence.random_sampling(
                    "ITERATION.fa"             ,
                    "SUBSAMPLE/SUBSAMPLE.fa"   ,
                    SUBSAMPLE_SIZE             ,
                    seed=SEED                  
                )
                shutil.rmtree ("EPANG")
                os    .mkdir  ("EPANG")
                shutil.rmtree ("TREE")
                os    .mkdir  ("TREE")
            else:
                rename_sequence.random_sampling(
                    INPUT_FA                   ,
                    "SUBSAMPLE/SUBSAMPLE.fa"   ,
                    SUBSAMPLE_SIZE             ,
                    seed=SEED                  ,
                    n = seq_count
                )
            
            #################
            #rename sequence#
            #################
            rename_sequence.rename_sequence(
                "SUBSAMPLE/SUBSAMPLE.fa"           , 
                "SUBSAMPLE/RENAMED_"+str(i)+".fa"
            )
            
            
            #######################################
            #construct subsample tree as reference#
            #######################################
            os.chdir(WD+"/TREE")
            subprocess.call(
                "bash " + CODEDIR+ "/shell/TREE.sh"             + 
                " -n "  + str(tree_thread_num)                  +
                " -m "  + TREEMETHOD                            +
                " -a "  + ALIGNED                               +
                " -f "  + WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa" +
                " -c "  + CODEDIR                               +
                " -w "  + WD+"/TREE"                            +
                " -p \""+ str(OPTION) + "\""                    +
                " -d "  + MODEL                                 + 
                " -q "  + SOFTWARE                              +
                " -s "  +ALIGNER                                ,
                shell=True
            )
            
            os.chdir(WD+"/PARAM")
            
            #if(seq_count>=30):
            if(raxml_thread_num>1):
                subprocess.call(
                    RAXMLPAR                                                        +
                    " -T "   + str(raxml_thread_num)                                +
                    " -f e"                                                         +
                    " -s "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned"      +
                    " -t "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned.tree" +
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            else:
                subprocess.call(
                    RAXMLSEQ                                                        +
                    " -f e"                                                         +
                    " -s "   + WD +"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned"       + 
                    " -t "   + WD +"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned.tree"  +
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            
            
            ########################################
            #Phylogenetic placement & visualization#
            ########################################
            os.chdir(WD)
            nodenum = (NODE_COUNT*seq_count)//INIT_SEQ_COUNT-1

            # select sequence file to place
            if(os.path.isfile(WD+"/INPUT.fa.aligned")):
                QUERY_FA = WD+"/INPUT.fa.aligned"
                ALIGNED_FOR_PLACEMENT = "aligned"
            else:
                QUERY_FA = WD+"/INPUT.fa"
                ALIGNED_FOR_PLACEMENT = ALIGNED

            # conduct placement
            placement.distributed_placement(
                WD                                                  , 
                EPANG                                               ,  
                WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned"       , 
                WD+"/PARAM/RAxML_result.PARAM_"+str(i)              , 
                WD+"/PARAM/RAxML_info.PARAM_"+str(i)                , 
                QUERY_FA                                            , 
                WD+"/EPANG"                                         , 
                THREAD_NUM                                          , 
                nodenum                                             ,
                CODEDIR                                             ,
                seq_count                                           ,
                ML_or_MP                                            ,
                RAXMLSEQ                                            ,
                ALIGNED_FOR_PLACEMENT                               ,
                SEED                                                ,
                hmm_aligner=HMM_ALIGNER                             ,
                hmm_profiler=HMM_PROFILER
            )
            
            
            ####################
            #parse .jplace file#
            ####################
            os.chdir(WD)
            para, Nseq_in_largest_subclade = \
                partition.partition(
                    WD+"/EPANG/placement_tree.out"          ,
                    WD+"/EPANG/edge_to_seqname_all.out"     ,
                    WD+"/PARTITION/partition"+str(i)+".out" ,
                    depth
                )
            print("detected "+str(para)+" paraphyletic sequences")
            
            
            ##################################################
            #get paraphyletic sequences and make ITERATION.fa#
            ##################################################
            if(para>prev_para or not(Nseq_in_largest_subclade<seq_count-1)):
                # The previous subsample achieved minimum number of paraphyletic groups
                i    -= 1 
                para  = prev_para
                break
            if(para!=0):
                
                # select subsample sequence file
                if os.path.isfile("EPANG/ref_query.fa.ref"):
                    ALIGNED_SUBSAMPLE = "EPANG/ref_query.fa.ref"
                else:
                    ALIGNED_SUBSAMPLE = "SUBSAMPLE/RENAMED_"+str(i)+".fa.aligned"

                shutil.copyfile(
                    ALIGNED_SUBSAMPLE,
                    "ITERATION.fa"
                )
                
                # select all sequence file
                if os.path.isfile(WD+"/INPUT.fa.aligned"):
                    ALIGNED_ALL = WD+"/INPUT.fa.aligned"
                    ALIGNED     = "aligned"
                else:
                    ALIGNED_ALL = "INPUT.fa" 

                partition.add_paraphyletic_fa(
                    WD+"/PARTITION/partition"+str(i)+".out" ,
                    "ITERATION.fa"                          ,
                    ALIGNED_ALL                             ,
                    SUBSAMPLE_SIZE                          ,
                    para
                    )
                i+=1
                prev_para=para
            else:
                break
        ##################
        #partition .fasta#
        ##################
        os.chdir(WD)
        print(
            "<Sequence count> Input:" + str(seq_count)                + "\n" +
            "  Largest subclade    :" + str(Nseq_in_largest_subclade) + "\n" +
            "  Problematic         :" + str(para)
        )
        
        if(i==-1): 
            print("Error: FRACluster.py cannot divide sequences into multiple subclades")
            sys.exit()
        
        FASTA_LIST=[WD+"/INPUT.fa"]
        if os.path.isfile(WD+"/INPUT.fa.aligned"): FASTA_LIST.append(WD+"/INPUT.fa.aligned")
        DIRdict = partition.partition_fasta(
            FASTA_LIST,
            NUMFILE,
            NODESDIR,
            WD,
            WD+"/PARTITION/partition"+str(min(i,MAX_ITERATION-1))+".out",
            "PARTITION.info",
            "UPSTREAM.nwk",
            WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa",
            ROOTING
            )
        
        partition.qsub_prep(
            ARGVS,
            QSUBDIR, 
            DIRdict, 
            INIT_SEQ_COUNT,
            seq_count_when_aligned
        )
        ##################
        #delete files    #
        ##################
        os.remove("INPUT.fa")
        os.remove("INPUT.aligned.fa")
        os.remove("tmp.nwk")
        '''
        shutil.rmtree("ANCSEQ")
        shutil.rmtree("EPANG")
        shutil.rmtree("PARAM")
        shutil.rmtree("PARTITION")
        shutil.rmtree("SUBSAMPLE")
        shutil.rmtree("TREE")
        '''
    
    elapsed_time=time.time()-start
    with open(WD+"/time.out", 'w') as handle:
        handle.write(
            str(seq_count)    + "," + 
            str(elapsed_time) + "," +
            "sec"             + "," +
            str(i)            + "," +
            "subsamplings"    + "\n"
        )

if __name__ == "__main__":
    argvs = sys.argv
    if (len(argvs)==23):
        FRACluster(
            argvs,
            argvs[1],
            int(argvs[2]),
            int(argvs[3]),
            argvs[4],
            int(argvs[5]),
            int(argvs[6]),
            argvs[7],
            argvs[8],
            argvs[9],
            argvs[10],
            argvs[11],
            argvs[12],
            argvs[13],
            argvs[14],
            argvs[15],
            argvs[16],
            argvs[17],
            argvs[18],
            int(argvs[19]),
            int(argvs[20]),
            argvs[21],
            argvs[22])
    elif ((len(argvs)==23+4)):
        FRACluster(
            argvs,
            argvs[1],
            int(argvs[2]),
            int(argvs[3]),
            argvs[4],
            int(argvs[5]),
            int(argvs[6]),
            argvs[7],
            argvs[8],
            argvs[9],
            argvs[10],
            argvs[11],
            argvs[12],
            argvs[13],
            argvs[14],
            argvs[15],
            argvs[16],
            argvs[17],
            argvs[18],
            int(argvs[19]),
            int(argvs[20]),
            argvs[21],
            argvs[22], 
            ALIGNER=argvs[23],
            HMM_PROFILER=argvs[24],
            HMM_ALIGNER=argvs[25],
            seq_count_when_aligned=int(argvs[26])
        )
    else:
        print("Error: Number of arguments: "+str(len(argvs))+" for FRACluster.py is wrong!")