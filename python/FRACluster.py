'''
While there are any problematic sequences,
    i) repeat 1-4 up to 5 times
        1. Random sampling (in : FASTA, out : FASTA)
        2. Phylogeny inference (in : FASTA, out : tree)
        3. RAxML parameter optimization (in : FASTA, tree, out RAxML log file)
        4. Phylogenetic placement (in : FASTA, Tree, RAxML log file, Query FASTA, out : .jplace file)
    ii) try to classify sequences into subclade (in : .jplace files, out : divided FASTAs, directory info file)
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
    if rename_sequence.outgroup_check_fast(WD+"/INPUT.fa.gz"):
        print("INPUT.fa properly include outgroup sequence")
    else:
        print("No sequence named \"root\" in INPUT.fa!\
               INPUT.fa should include outgroup sequence named \"root\"")
        return
    
    # check number of sequences
    seq_array                 = rename_sequence.count_sequence_fast("INPUT.fa.gz")
    seq_count                 = seq_array[0]
    if(INIT_SEQ_COUNT==0): 
        INIT_SEQ_COUNT        = seq_count # only in d0
        seq_count_when_aligned= None
    seq_length                = seq_array[1]
    raxml_thread_num          = min(max(seq_length//500,2),THREAD_NUM) # use 1 thread per 500 bp in RAxML
    depth                     = max(math.floor(math.log2(seq_count/THRESHOLD))+2,2)
    tree_thread_num           = THREAD_NUM

    # check if aligned
    if (ALIGNED=="unaligned"):
        if(os.path.isfile(WD+"/INPUT.fa.aligned.gz")):
            if (seq_count < seq_count_when_aligned * ALIGNMENT_TIMING_PARAMETER):

                print(WD,seq_count, seq_count_when_aligned, "alignment needed!!")

                os.remove(WD+"/INPUT.fa.aligned.gz")
                INPUT_FA               = WD+"/INPUT.fa.gz"
                seq_count_when_aligned = seq_count
                
            else:
                INPUT_FA = WD+"/INPUT.fa.aligned.gz"
                ALIGNED  = "aligned"
        else:

            print(WD,seq_count,seq_count_when_aligned,"alignment needed!!")

            INPUT_FA               = WD+"/INPUT.fa.gz"
            seq_count_when_aligned = seq_count
    elif (ALIGNED=="aligned"):
        INPUT_FA               = WD+"/INPUT.fa.gz"

    i=1
    
    # call direct tree reconstruction

    if(seq_count<=THRESHOLD):
        if(seq_count<4):
            partition.tiny_tree(INPUT_FA,"TERMINAL.nwk")
            print("seq_count < 4!")
        else:
            os.mkdir(    "TREE")
            os.chdir(WD+"/TREE")
            subprocess.call("bash " + CODEDIR + "/shell/TREE.sh" +
                            " -n "  + str(tree_thread_num)       +
                            " -m "  + TREEMETHOD                 +
                            " -a "  + ALIGNED                    +
                            " -f "  + INPUT_FA                   +
                            " -c "  + CODEDIR                    +
                            " -w "  + WD+"/TREE"                 +
                            " -p \""+ str(OPTION)+"\""           +
                            " -d "  + MODEL                      +
                            " -q "  + SOFTWARE                   +
                            " -s "  + ALIGNER                    ,
                            shell   = True                       )
            partition.rooting_and_remove(
                INPUT_FA+".aligned.tree"     ,
                WD+"/TERMINAL.nwk"           ,
                "root"
            )
    
    # call fractal FRACTAL # don't forget to change!!!!!
    elif(NODE_COUNT>1 and seq_count<=INIT_SEQ_COUNT//NODE_COUNT): 
        # Quit distribution after the available computer node saturated 
        os.mkdir("TREE")
        os.chdir(WD+"/TREE")
        FRACTAL_COMMAND = "FRACTAL"                         + \
                          " -i " + INPUT_FA                 + \
                          " -k " + str(SUBSAMPLE_SIZE)      + \
                          " -b " + MODEL                    + \
                          " -p " + ML_or_MP                 + \
                          " -t " + str(THRESHOLD)           + \
                          " -x " + str(MAX_ITERATION)       + \
                          " -c " + str(THREAD_NUM)          + \
                          " -r " + SEED                     + \
                          " -e "
        
        if (TREEMETHOD!="unspecified"): 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -m "+TREEMETHOD
        else: 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -s "+SOFTWARE
        if (ALIGNED=='unaligned'):
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -u "
        if (OPTION!=""):
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -a " + OPTION
        
        subprocess.call(FRACTAL_COMMAND,shell=True)
        shutil.move(WD+"/TREE/FRACTALout.nwk",WD+"/TERMINAL.nwk")

    # call FRACTAL cycle
    else:
        os.mkdir("SUBSAMPLE")
        os.mkdir("TREE")
        os.mkdir("PARAM")
        os.mkdir("EPANG")
        os.mkdir("PARTITION")
        
        prev_para = seq_count

        while i<MAX_ITERATION:
            
            os.chdir(WD)

            #################
            #random sampling#
            #################
            if(os.path.isfile(WD+"/INPUT.fa.aligned.gz")):
                INPUT_FA = WD+"/INPUT.fa.aligned"
                ALIGNED  = "aligned" 

            if(os.path.isfile("ITERATION.fa.gz.gz")):
                rename_sequence.random_sampling(
                    "ITERATION.fa.gz"             ,
                    "SUBSAMPLE/SUBSAMPLE.fa.gz"   ,
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
                    "SUBSAMPLE/SUBSAMPLE.fa.gz"   ,
                    SUBSAMPLE_SIZE             ,
                    seed=SEED                  ,
                    n = seq_count
                )
            
            #################
            #rename sequence#
            #################
            rename_sequence.rename_sequence(
                "SUBSAMPLE/SUBSAMPLE.fa.gz"           , 
                "SUBSAMPLE/RENAMED_"+str(i)+".fa.gz"
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
                " -f "  + WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz" +
                " -c "  + CODEDIR                               +
                " -w "  + WD+"/TREE"                            +
                " -p \""+ str(OPTION) + "\""                    +
                " -d "  + MODEL                                 + 
                " -q "  + SOFTWARE                              +
                " -s "  + ALIGNER                               ,
                shell=True
            )
            
            os.chdir(WD+"/PARAM")
            
            if(raxml_thread_num>1):
                subprocess.call(
                    RAXMLPAR                                                        +
                    " -T "   + str(raxml_thread_num)                                +
                    " -f e"                                                         +
                    " -s "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"      +
                    " -t "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree" +
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            else:
                subprocess.call(
                    RAXMLSEQ                                                        +
                    " -f e"                                                         +
                    " -s "   + WD +"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"       + 
                    " -t "   + WD +"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree"  +
                    " -n "   + "PARAM_"+str(i)                                      +
                    # " -n "   + "PARAM_"+str(i-i%2)                                  + # for test
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            
            # if parameter optimization failed
            if (not os.path.isfile(WD+"/PARAM/RAxML_result.PARAM_"+str(i))):

                print("Parameter optimization failed...")

                # if i == 0, start from random sampling again, else use the result of previous i
                if(i > 1 and os.path.isfile(WD+"/PARTITION/partition"+str(i-1)+".out") ): 
                    i -= 1
                    para  = prev_para
                    break
                else:
                    if(os.path.isfile(WD+"/ITERATION.fa.gz")):
                        os.remove(WD+"/ITERATION.fa.gz")
                    i += 1
            
            else: # if parameter optimization succeeded

                ########################################
                #Phylogenetic placement & visualization#
                ########################################
                os.chdir(WD)
                nodenum = (NODE_COUNT*seq_count)//INIT_SEQ_COUNT-1

                # select sequence file to place
                if(os.path.isfile(WD+"/INPUT.fa.gz.aligned")):
                    QUERY_FA = WD+"/INPUT.fa.gz.aligned"
                    ALIGNED_FOR_PLACEMENT = "aligned"
                else:
                    QUERY_FA = WD+"/INPUT.fa"
                    ALIGNED_FOR_PLACEMENT = ALIGNED

                # conduct placement
                placement.distributed_placement(
                    WD                                                  , 
                    EPANG                                               ,  
                    WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"       , 
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
                #get paraphyletic sequences and make ITERATION.fa.gz#
                ##################################################
                if(para>prev_para or not(Nseq_in_largest_subclade<seq_count-1)):
                    
                    # The previous subsample achieved minimum number of paraphyletic groups
                    i    -= 1 
                    para  = prev_para
                    break

                if(para!=0): # if problematic sequences remained

                    # select subsample sequence file
                    if os.path.isfile(WD+"/SUBSAMPLE.fa.aligned"):
                        ALIGNED_SUBSAMPLE = WD+"/SUBSAMPLE.fa.gz.aligned"
                    else:
                        ALIGNED_SUBSAMPLE = "SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"

                    shutil.copyfile(
                        ALIGNED_SUBSAMPLE,
                        "ITERATION.fa.gz"
                    )
                    
                    # select all sequence file
                    if os.path.isfile(WD+"/INPUT.fa.aligned"):
                        ALIGNED_ALL = WD+"/INPUT.fa.aligned"
                        ALIGNED     = "aligned"
                    else:
                        ALIGNED_ALL = "INPUT.fa" 

                    partition.add_paraphyletic_fa(
                        WD+"/PARTITION/partition"+str(i)+".out" ,
                        "ITERATION.fa.gz"                          ,
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
        try:
            print(
                "<Sequence count> Input:" + str(seq_count)                + "\n" +
                "  Largest subclade    :" + str(Nseq_in_largest_subclade) + "\n" +
                "  Problematic         :" + str(para)
            )
        except:
            None
        
        if(i==-1): 
            print("Error: FRACluster.py cannot divide sequences into multiple subclades")
            sys.exit()
        
        FASTA_LIST = [ WD+"/INPUT.fa" ]
        if os.path.isfile(WD+"/INPUT.fa.gz.aligned"):
            FASTA_LIST.append(WD+"/INPUT.fa.gz.aligned")
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

        try: 
            shutil.rmtree("TREE")
        except:
            None

    ##################
    #delete files    #
    ##################

    os.chdir(WD)

    filenames = [
        "INPUT.fa",
        "ITERATION.fa.gz",
        "tmp.nwk",
        "INPUT.fa.aligned",
        "SUBSAMPLE.fa.aligned"
        ]
    
    dirnames = [
        "EPANG",
        "PARAM",
        "PARTITION",
        "SUBSAMPLE",
        #"TREE"
        ]

    for filename in filenames:
        try:
            os.remove(filename)
        except:
            None

    for dirname in dirnames:
        try: 
            shutil.rmtree(dirname)
        except:
            None
        
    
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

    if (len(argvs)==24):
        FRACluster(
            ARGVS             = argvs, 
            WD                = argvs[1],
            MAX_ITERATION     = int(argvs[2]),
            SUBSAMPLE_SIZE    = int(argvs[3]),   
            NODESDIR          = argvs[4],
            THRESHOLD         = int(argvs[5]), 
            THREAD_NUM        = int(argvs[6]),
            NUMFILE           = argvs[7],
            QSUBDIR           = argvs[8], 
            CODEDIR           = argvs[9],
            ROOTING           = argvs[10],
            MODEL             = argvs[11],
            OPTION            = argvs[12],
            TREEMETHOD        = argvs[13],
            ALIGNED           = argvs[14],
            EPANG             = argvs[15],
            RAXMLSEQ          = argvs[16],
            RAXMLPAR          = argvs[17],
            SOFTWARE          = argvs[18],
            NODE_COUNT        = int(argvs[19]),
            INIT_SEQ_COUNT    = int(argvs[20]),
            SEED              = argvs[21],
            ML_or_MP          = argvs[22],
            ALIGNER           = "unspecified", 
            HMM_PROFILER      = "unspecified", 
            HMM_ALIGNER       = "unspecified",
            seq_count_when_aligned=None
            )
    elif ((len(argvs)==23+4)):
        FRACluster(
            ARGVS             = argvs, 
            WD                = argvs[1],
            MAX_ITERATION     = int(argvs[2]),
            SUBSAMPLE_SIZE    = int(argvs[3]),   
            NODESDIR          = argvs[4],
            THRESHOLD         = int(argvs[5]), 
            THREAD_NUM        = int(argvs[6]),
            NUMFILE           = argvs[7],
            QSUBDIR           = argvs[8], 
            CODEDIR           = argvs[9],
            ROOTING           = argvs[10],
            MODEL             = argvs[11],
            OPTION            = argvs[12],
            TREEMETHOD        = argvs[13],
            ALIGNED           = argvs[14],
            EPANG             = argvs[15],
            RAXMLSEQ          = argvs[16],
            RAXMLPAR          = argvs[17],
            SOFTWARE          = argvs[18],
            NODE_COUNT        = int(argvs[19]),
            INIT_SEQ_COUNT    = int(argvs[20]),
            SEED              = argvs[21],
            ML_or_MP          = argvs[22],
            ALIGNER           = argvs[23],
            HMM_PROFILER      = argvs[24],
            HMM_ALIGNER       = argvs[25],
            seq_count_when_aligned=int(argvs[26])
        )
    else:
        print("Error: Number of arguments: "+str(len(argvs))+" for FRACluster.py is wrong!")