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
import extraction
import error_process

import math
import time
import random

def FRACluster(ARGVS, WD, MAX_ITERATION, SUBSAMPLE_SIZE, NODESDIR, THRESHOLD, THREAD_NUM, NUMFILE, QSUBDIR, CODEDIR, 
               ROOTING, MODEL, OPTION,TREEMETHOD, ALIGNED, EPANG, RAXMLSEQ, RAXMLPAR, SOFTWARE,NODE_COUNT,
               INIT_SEQ_COUNT,SEED,ML_or_MP, EXTRACTION_SIZE,careful,
               ALIGNER="unspecified", HMM_PROFILER="unspecified", HMM_ALIGNER="unspecified",
               seq_count_when_aligned=None,
               ):
    
    # start timer
    start = time.time() 

    ##### hyper parameter #####
    ALIGNMENT_TIMING_PARAMETER = 0.5
    SPLIT_THRESHOLD            = 10000
    ###########################

    ### get input file name ###
    infile_name = os.listdir(WD)[0]
    infile_path = WD+"/"+infile_name
    infile_aligned_path = WD+"/"+".".join(infile_name.split(".")[:-1]+["aligned",infile_name.split(".")[-1]])
    ###########################
    
    ## check input file property ##
    try:
        root_idx = rename_sequence.outgroup_check_fast(infile_path, "fasta")
    except:
        error_process.no_root(); exit(1)
    seq_count                 = rename_sequence.count_sequence_fast(infile_path)
    is_gzipped                = (infile_path.split(".")[-1] == "gz")
    if (is_gzipped):
        gzip_extention        = ".gz"
    else:
        gzip_extention        = ""
    if(INIT_SEQ_COUNT==0): 
        INIT_SEQ_COUNT        = seq_count # only in d0
        seq_count_when_aligned= None
    # check if aligned
    if (ALIGNED=="unaligned"):

        if(os.path.isfile(infile_aligned_path)):
            if (seq_count < seq_count_when_aligned * ALIGNMENT_TIMING_PARAMETER):

                print(WD,seq_count, seq_count_when_aligned, "alignment needed!!")

                os.remove(infile_aligned_path)
                INPUT_FA               = infile_path
                seq_count_when_aligned = seq_count
                
            else:
                INPUT_FA = infile_aligned_path
                ALIGNED  = "aligned"
        else:
            print(
                "Directory:",WD,
                "#sequence:",seq_count,
                "#sequence in the last alignment:",seq_count_when_aligned,
                "alignment needed!!"
                )
            INPUT_FA               = infile_path
            seq_count_when_aligned = seq_count
    elif (ALIGNED=="aligned"):
        INPUT_FA               = infile_path
    ################################

    ######## parameter ########
    if(SEED=="random"): random.seed(int(random.randint(0,99999)))
    elif(len(SEED)!=0): random.seed(int(SEED))
    else:print("-r Error: invalid random seed!")
    raxml_thread_num          = 1
    depth                     = max(math.floor(math.log2(seq_count/THRESHOLD))+2,2)
    tree_thread_num           = THREAD_NUM
    ###########################

    # move to Working Directory
    os.chdir(WD) 

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
        if (is_gzipped):
            gzip_option = " -g "
        else:
            gzip_option = ""
        FRACTAL_COMMAND = "FRACTAL"                         + \
                          " -i " + INPUT_FA                 + \
                          " -k " + str(SUBSAMPLE_SIZE)      + \
                          " -b " + MODEL                    + \
                          " -p " + ML_or_MP                 + \
                          " -t " + str(THRESHOLD)           + \
                          " -x " + str(MAX_ITERATION)       + \
                          " -c " + str(THREAD_NUM)          + \
                          " -r " + SEED                     + \
                          " -e "                            + \
                          gzip_option
        
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
        nodenum = (NODE_COUNT*seq_count)//INIT_SEQ_COUNT-1


        while i<MAX_ITERATION:
            
            os.chdir(WD)

            ### get input file name ###
            iterationfile_path         = WD + "ITERATION.fa"                    + gzip_extention
            subsamplefile_path         = WD + "SUBSAMPLE/SUBSAMPLE.fa"          + gzip_extention
            renamed_subsamplefile_path = WD + "SUBSAMPLE/RENAMED_"+str(i)+".fa" + gzip_extention
            ###########################

            ################
            # split or not #
            ################
            split = seq_count > SPLIT_THRESHOLD or nodenum > 1

            #################
            #random sampling#
            #################
            if(os.path.isfile(infile_aligned_path)):
                INPUT_FA = infile_aligned_path
                ALIGNED  = "aligned"

            if (i>1):
                shutil.rmtree ("EPANG")
                os    .mkdir  ("EPANG")
                shutil.rmtree ("TREE")
                os    .mkdir  ("TREE")

            if(os.path.isfile(iterationfile_path)):
                if (split): # if split
                    if os.path.isfile(infile_aligned_path) and not os.path.exists(infile_aligned_path+".split"):
                        Nseq_per_file = min(10000, seq_count//max(nodenum,1))
                        subprocess.call("seqkit split2 -s "+str(Nseq_per_file)+" "+infile_aligned_path+" &> /dev/null", shell=True)
                sampled_seq_name_list = \
                    rename_sequence.random_sampling(
                        iterationfile_path            ,
                        subsamplefile_path            ,
                        SUBSAMPLE_SIZE                ,
                        seed=SEED
                    )
            else:
                # decompose FASTA
                if (seq_count > 10000 or nodenum > 1):
                    Nseq_per_file = min(10000, seq_count//max(nodenum,1))
                    subprocess.call("seqkit split2 -s "+str(Nseq_per_file)+" "+INPUT_FA, shell=True)
                    # subsampling
                    sampled_seq_name_list = \
                        rename_sequence.random_sampling_from_splitted(
                            INPUT_FA+".split"          ,
                            subsamplefile_path         ,
                            SUBSAMPLE_SIZE             ,
                            SEED                       ,
                            Nseq_per_file              ,
                            root_idx                   ,
                            n = seq_count
                        )
                else:
                    # subsampling
                    sampled_seq_name_list = \
                        rename_sequence.random_sampling(
                            INPUT_FA                   ,
                            subsamplefile_path         ,
                            SUBSAMPLE_SIZE             ,
                            seed=SEED                  ,
                            n = seq_count
                        )
            
            #################
            #rename sequence#
            #################
            seqname2renamedname = rename_sequence.rename_sequence(
                subsamplefile_path, 
                renamed_subsamplefile_path
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
                " -f "  + renamed_subsamplefile_path            +
                " -c "  + CODEDIR                               +
                " -w "  + WD+"/TREE"                            +
                " -p \""+ str(OPTION) + "\""                    +
                " -d "  + MODEL                                 + 
                " -q "  + SOFTWARE                              +
                " -s "  + ALIGNER                               ,
                shell=True
            )
            
            # Subsample tree extraction

            if ( EXTRACTION_SIZE < len(sampled_seq_name_list)):
                
                extracted_seq_name_list   = ["root"]+random.sample(sampled_seq_name_list, EXTRACTION_SIZE)
                extracted_seq_rename_list = [ seqname2renamedname[name] for name in extracted_seq_name_list ]
                extraction.tree_extraction(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree",
                    set(extracted_seq_rename_list)                         ,
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.extracted.tree"
                    )
                extraction.fasta_extraction(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz"          ,
                    set(extracted_seq_rename_list)                      ,
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz",
                    )
                subprocess.call(
                    "cat "                                                          +
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz | "         +
                    "gunzip > " + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted" ,
                    shell=True
                )
                os.remove(WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree")
                os.remove(WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz")
                os.remove(WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned")
                
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.extracted.tree",
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree"
                    )
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz",
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz"
                )
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted",
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"
                )

            os.chdir(WD+"/PARAM")
            
            if(raxml_thread_num>1):
                subprocess.call(
                    RAXMLPAR                                                        +
                    " -T "   + str(raxml_thread_num)                                +
                    " -f e"                                                         +
                    " -s "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"   +
                    " -t "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree"+
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            else:
                subprocess.call(
                    RAXMLSEQ                                                        +
                    " -f e"                                                         +
                    " -s "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"   + 
                    " -t "   + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned.tree"+
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

                # select sequence file to place
                if(os.path.isfile(WD+"/INPUT.fa.gz.aligned")):
                    QUERY_FA = WD+"/INPUT.fa.gz.aligned"
                    ALIGNED_FOR_PLACEMENT = "aligned"
                else:
                    QUERY_FA = WD+"/INPUT.fa.gz"
                    ALIGNED_FOR_PLACEMENT = ALIGNED
                
                if (ALIGNED=="aligned"):
                    refseq = WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz"
                else:
                    refseq = WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"

                # conduct placement
                placement.distributed_placement(
                    WD                                                  , 
                    EPANG                                               ,  
                    refseq                                              , 
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
                    careful=careful                                     ,
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

                #####################################################
                #get paraphyletic sequences and make ITERATION.fa.gz#
                #####################################################
                resampling_needed = False
                if(para>prev_para):
                    
                    # if i == 0, start from random sampling again, else use the result of previous i
                    if(i > 1 and os.path.isfile(WD+"/PARTITION/partition"+str(i-1)+".out") ): 
                        i -= 1
                        para  = prev_para
                        break
                    else:
                        if(os.path.isfile(WD+"/ITERATION.fa.gz")):
                            os.remove(WD+"/ITERATION.fa.gz")
                        i += 1
                        resampling_needed = True
                
                if( not Nseq_in_largest_subclade < seq_count-1 ):
                    
                    i += 1
                    resampling_needed = True

                if(para!=0): # if problematic sequences remained

                    # select subsample sequence file
                    if os.path.isfile(WD+"/SUBSAMPLE.fa.aligned"):
                        ALIGNED_SUBSAMPLE = WD+"/SUBSAMPLE.fa.aligned"
                    else:
                        ALIGNED_SUBSAMPLE = "SUBSAMPLE/RENAMED_"+str(i)+".fa.gz.aligned"

                    #shutil.copyfile(
                    #    ALIGNED_SUBSAMPLE,
                    #    "ITERATION.fa.gz"
                    #)
                    subprocess.call(
                        "gzip -c "+ALIGNED_SUBSAMPLE+" > "+WD+"/ITERATION.fa.gz",
                        shell=True
                    )
                    
                    # select all sequence file
                    if os.path.isfile(WD+"/INPUT.fa.gz.aligned"):
                        ALIGNED_ALL = WD+"/INPUT.fa.gz.aligned"
                        ALIGNED     = "aligned"
                    else:
                        ALIGNED_ALL = "INPUT.fa.gz" 

                    partition.add_paraphyletic_fa(
                        WD+"/PARTITION/partition"+str(i)+".out" ,
                        "ITERATION.fa.gz"                       ,
                        ALIGNED_ALL                             ,
                        SUBSAMPLE_SIZE                          ,
                        para
                        )
                    i+=1
                    prev_para=para
                elif (not resampling_needed):
                    print("resampling",WD,Nseq_in_largest_subclade, seq_count, resampling_needed,i,"break")
                    break
                else:
                    prev_para = 0
                    print("resampling",WD,Nseq_in_largest_subclade, seq_count, resampling_needed, i)

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
        
        if (Nseq_in_largest_subclade == seq_count - 1 ): # if all sequences were classified into one subclade, FRACTAL gives up for inference of this clade
            print ("Error: FRACluster.py could not divide sequences into multiple subclades")
            return
        
        FASTA_LIST = [ WD+"/INPUT.fa.gz" ]
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
            WD+"/SUBSAMPLE/RENAMED_"+str(i)+".fa.gz",
            ROOTING,
            nodenum = nodenum,
            codedir = CODEDIR
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
        "INPUT.fa.gz",
        "INPUT.fa.aligned",
        "INPUT.fa.aligned.gz",
        "ITERATION.fa.gz",
        "tmp.nwk",
        "SUBSAMPLE.fa.aligned.gz",
        "INPUT.fa.gz.aligned",
        "INPUT.fa.gz.aligned.tree"
        ]
    
    dirnames = [
        "INPUT.fa.gz.split",
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

    print(argvs)

    if (len(argvs)==26):
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
            EXTRACTION_SIZE   = int(argvs[23]),  
            careful           = int(argvs[24]),
            ALIGNER           = "unspecified", 
            HMM_PROFILER      = "unspecified", 
            HMM_ALIGNER       = "unspecified",
            seq_count_when_aligned=None
            )
    elif ((len(argvs)==26+3)):
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
            EXTRACTION_SIZE   = int(argvs[23]),
            careful           = int(argvs[24]),
            ALIGNER           = argvs[25],
            HMM_PROFILER      = argvs[26],
            HMM_ALIGNER       = argvs[27],
            seq_count_when_aligned=int(argvs[28])
        )
    else:
        print("Error: Number of arguments: "+str(len(argvs))+" for FRACluster.py is wrong!")