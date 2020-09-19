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

    # remove empty input files#
    for infile_name in os.listdir(WD):
        if (os.path.getsize(WD+"/"+infile_name) == 0):
            os.remove(WD+"/"+infile_name)
    ###########################

    ### get input file name ###
    infile_namelist              = os.listdir(WD)
    infile_pathlist              = [WD+"/"+name for name in infile_namelist]
    example_infile_fpath         = infile_pathlist[0]
    infile_pathlist_aligned      = [infile_path+".aligned" for infile_path in infile_pathlist]
    example_infile_fpath_aligned = infile_pathlist_aligned[0]
    iterationfile_path="unspecified"
    ###########################
    
    ## check input file property ##
    try:
        root_fpath = rename_sequence.outgroup_check_fast(infile_pathlist, "fasta")
    except:
        error_process.no_root(); exit(1)
    fpath2seqcount            = rename_sequence.count_sequence_fast(infile_pathlist)
    seq_count                 = sum(fpath2seqcount.values())
    is_gzipped                = (example_infile_fpath.split(".")[-1] == "gz")
    if (is_gzipped):
        gzip_extention        = ".gz"
        gzip_command          = "gzip"
    else:
        gzip_extention        = ""
        gzip_command          = "cat"
    if(INIT_SEQ_COUNT==0): 
        INIT_SEQ_COUNT        = seq_count # only in d0
        seq_count_when_aligned= None
    # check if aligned
    '''
    if (ALIGNED=="unaligned"):

        if(os.path.isfile(infile_aligned_path)):
            if (seq_count < seq_count_when_aligned * ALIGNMENT_TIMING_PARAMETER):

                print(WD,seq_count, seq_count_when_aligned, "alignment needed!!")

                os.remove(infile_aligned_path)
                seq_count_when_aligned = seq_count
                
            else:
                ALIGNED  = "aligned"
        else:
            print(
                "Directory:",WD,
                "#sequence:",seq_count,
                "#sequence in the last alignment:",seq_count_when_aligned,
                "alignment needed!!"
                )
            seq_count_when_aligned = seq_count
    elif (ALIGNED=="aligned"):
        INPUT_FA               = infile_path
    '''
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
        if len(infile_pathlist)>1:
            concat_infpath = WD+"/INPUT.fa"+gzip_extention
            subprocess.call(
                "cat "+" ".join(infile_pathlist)+"|"+gzip_command+" > "+concat_infpath,
                shell=True
            )
        else:
            concat_infpath = infile_pathlist[0]
        if(seq_count<4):
            partition.tiny_tree(concat_infpath,"TERMINAL.nwk")
            print("seq_count < 4!")
        else:
            os.mkdir(    "TREE")
            os.chdir(WD+"/TREE")
            subprocess.call("bash " + CODEDIR + "/shell/TREE.sh" +
                            " -n "  + str(tree_thread_num)       +
                            " -m "  + TREEMETHOD                 +
                            " -a "  + ALIGNED                    +
                            " -f "  + concat_infpath             +
                            " -c "  + CODEDIR                    +
                            " -w "  + WD+"/TREE"                 +
                            " -p \""+ str(OPTION)+"\""           +
                            " -d "  + MODEL                      +
                            " -q "  + SOFTWARE                   +
                            " -s "  + ALIGNER                    ,
                            shell   = True                       )
            partition.rooting_and_remove(
                concat_infpath+".aligned.tree" ,
                WD+"/TERMINAL.nwk"             ,
                "root"
            )
    
    # call fractal FRACTAL # don't forget to change!!!!!
    elif(NODE_COUNT>1 and seq_count<=INIT_SEQ_COUNT//NODE_COUNT):
        # To Do: remove concatenation
        concat_infpath = WD+"/INPUT.fa"+gzip_extention
        subprocess.call(
            "cat "+" ".join(infile_pathlist)+"|"+gzip_command+" > "+concat_infpath,
            shell=True
        ) 
        # Quit distribution after the available computer node saturated 
        os.mkdir("TREE")
        os.chdir(WD+"/TREE")
        if (is_gzipped):
            gzip_option = " -g "
        else:
            gzip_option = ""
        FRACTAL_COMMAND = "FRACTAL"                         + \
                          " -i " + concat_infpath           + \
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

            ### get input file name2###
            iterationfile_path         = WD + "/ITERATION.fa"                    + gzip_extention
            subsamplefile_path         = WD + "/SUBSAMPLE/SUBSAMPLE.fa"          
            renamed_subsamplefile_path = WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa" 
            ###########################
        
            ###########################
            # renew several directory #
            ###########################
            if (i>1):
                shutil.rmtree ("EPANG")
                os    .mkdir  ("EPANG")
                shutil.rmtree ("TREE")
                os    .mkdir  ("TREE")

            ############
            # set mode #   
            ############
            AFTER_ALIGNMENT = os.path.exists(example_infile_fpath_aligned+".split")
            if(AFTER_ALIGNMENT):
                ALIGNED         = "aligned"

            ################
            #  file split  #
            ################
            split = nodenum > 1
            Nseq_per_file = min(SPLIT_THRESHOLD, seq_count//max(nodenum,1))
            if (split):
                if (AFTER_ALIGNMENT):
                    file_pathlist_to_be_splitted = infile_pathlist_aligned
                else:
                    file_pathlist_to_be_splitted = infile_pathlist

                splitted_dirpath = file_pathlist_to_be_splitted[0]+".split"
                if not os.path.exists(splitted_dirpath):
                    for j, file_path in enumerate(file_pathlist_to_be_splitted):
                        subprocess.call("seqkit split2 -s "+str(Nseq_per_file)+" "+file_path+" &> /dev/null", shell=True)
                        splitted_fnamelist = sorted(os.listdir(file_path+".split"))
                        for k, splitted_fname in enumerate(splitted_fnamelist):
                            if (k < len(splitted_fnamelist)-1):
                                fpath2seqcount[file_pathlist_to_be_splitted[0]+".split/" + splitted_fname] = Nseq_per_file
                            else:
                                fpath2seqcount[file_pathlist_to_be_splitted[0]+".split/" + splitted_fname] = seq_count % Nseq_per_file
                    for file_path in file_pathlist_to_be_splitted[1:]:
                        subprocess.call(
                            "mv " + file_path+".split/* " + splitted_dirpath +
                            "rm -r " + file_path+".split",
                            shell=True
                        )
            else:
                os.mkdir(example_infile_fpath+".split")
                for infile_path in infile_pathlist:
                    fpath2seqcount[example_infile_fpath+".split/"+infile_path.split("/")[-1]] = fpath2seqcount[infile_path]
                    shutil.move(infile_path, example_infile_fpath+".split")
                splitted_dirpath = example_infile_fpath+".split"

            # rename
            for i, splitted_fname in enumerate(os.listdir(splitted_dirpath)):
                splitted_fpath   = splitted_dirpath+"/"+splitted_fname
                renamed_filepath = splitted_dirpath+"/INPUT.part"+str(i)+".fa"+gzip_extention
                shutil.move(splitted_fpath, renamed_filepath)
                fpath2seqcount[renamed_filepath] = fpath2seqcount[splitted_fpath]

            #################
            #random sampling#
            #################

            if (os.path.exists(iterationfile_path)):
                sampled_seq_name_list = \
                    rename_sequence.random_sampling_fasta( # fasta only
                        in_dirpath = None,
                        out_fname  = subsamplefile_path,
                        subsample_size = SUBSAMPLE_SIZE,
                        fpath2seqcount = None,
                        root_fpath     = root_fpath,
                        total_seqcount = None, 
                        file_format    = "fasta",
                        in_fpath       = iterationfile_path
                        )
            else:
                # subsampling
                sampled_seq_name_list = \
                    rename_sequence.random_sampling_fasta( # fasta only
                        in_dirpath = splitted_dirpath,
                        out_fname  = subsamplefile_path,
                        subsample_size = SUBSAMPLE_SIZE,
                        fpath2seqcount = fpath2seqcount,
                        root_fpath     = root_fpath,
                        total_seqcount = seq_count, 
                        file_format    = "fasta",
                        in_fpath       = None
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
                    renamed_subsamplefile_path+".aligned.tree",
                    set(extracted_seq_rename_list)                         ,
                    renamed_subsamplefile_path+".aligned.extracted.tree"
                    )
                extraction.fasta_extraction(
                    renamed_subsamplefile_path          ,
                    set(extracted_seq_rename_list)                      ,
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz",
                    )
                subprocess.call(
                    "cat "                                                          +
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz | "         +
                    "gunzip > " + WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted" ,
                    shell=True
                )
                os.remove(renamed_subsamplefile_path+".aligned.tree")
                os.remove(renamed_subsamplefile_path+"")
                os.remove(renamed_subsamplefile_path+".aligned")
                
                os.rename(
                    renamed_subsamplefile_path+".aligned.extracted.tree",
                    renamed_subsamplefile_path+".aligned.tree"
                    )
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.gz",
                    renamed_subsamplefile_path
                )
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted",
                    renamed_subsamplefile_path+".aligned"
                )

            os.chdir(WD+"/PARAM")
            
            if(raxml_thread_num>1):
                subprocess.call(
                    RAXMLPAR                                                        +
                    " -T "   + str(raxml_thread_num)                                +
                    " -f e"                                                         +
                    " -s "   + renamed_subsamplefile_path+".aligned"   +
                    " -t "   + renamed_subsamplefile_path+".aligned.tree"+
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                ,
                    shell=True
                )
            else:
                subprocess.call(
                    RAXMLSEQ                                                        +
                    " -f e"                                                         +
                    " -s "   + renamed_subsamplefile_path+".aligned"   + 
                    " -t "   + renamed_subsamplefile_path+".aligned.tree"+
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
                    if(os.path.isfile(iterationfile_path)):
                        os.remove(iterationfile_path)
                    i += 1
            
            else: # if parameter optimization succeeded

                ########################################
                #Phylogenetic placement & visualization#
                ########################################
                os.chdir(WD)

                # select sequence file to place
                if(os.path.isfile(example_infile_fpath_aligned+".split")):
                    QUERY_FA_DIR          = example_infile_fpath_aligned+".split"
                    ALIGNED_FOR_PLACEMENT = "aligned"
                else:
                    QUERY_FA_DIR          = example_infile_fpath+".split"
                    ALIGNED_FOR_PLACEMENT = ALIGNED
                
                if (ALIGNED=="aligned"):
                    refseq = renamed_subsamplefile_path
                else:
                    refseq = renamed_subsamplefile_path+".aligned"

                # conduct placement
                placement.distributed_placement(
                    WD                                                  , 
                    EPANG                                               ,  
                    refseq                                              , 
                    WD+"/PARAM/RAxML_result.PARAM_"+str(i)              , 
                    WD+"/PARAM/RAxML_info.PARAM_"+str(i)                , 
                    QUERY_FA_DIR                                        , 
                    WD+"/EPANG"                                         , 
                    THREAD_NUM                                          , 
                    nodenum                                             ,
                    CODEDIR                                             ,
                    seq_count                                           ,
                    ML_or_MP                                            ,
                    RAXMLSEQ                                            ,
                    ALIGNED_FOR_PLACEMENT                               ,
                    SEED                                                ,
                    alignment_outdir = example_infile_fpath_aligned+".split",
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
                        if(os.path.isfile(iterationfile_path)):
                            os.remove(iterationfile_path)
                        i += 1
                        resampling_needed = True
                
                if( not Nseq_in_largest_subclade < seq_count-1 ):
                    
                    i += 1
                    resampling_needed = True

                if(para!=0): # if problematic sequences remained

                    # select subsample sequence file
                    if os.path.isfile(subsamplefile_path+".aligned"):
                        ALIGNED_SUBSAMPLE = subsamplefile_path+".aligned"
                    else:
                        ALIGNED_SUBSAMPLE = subsamplefile_path
                        
                    subprocess.call(
                        "cat "+ALIGNED_SUBSAMPLE+" | "+gzip_command+"> "+iterationfile_path,
                        shell=True
                    )
                    
                    # select all sequence file
                    if os.path.isfile(infile_path+".aligned"):
                        ALIGNED_ALL = infile_path+".aligned"
                        ALIGNED     = "aligned"
                    else:
                        ALIGNED_ALL = infile_path+gzip_extention

                    partition.add_paraphyletic_fa(
                        WD+"/PARTITION/partition"+str(i)+".out" ,
                        iterationfile_path                      ,
                        ALIGNED_ALL                             ,
                        SUBSAMPLE_SIZE                          ,
                        para
                        )
                    i+=1
                    prev_para=para
                elif (not resampling_needed):
                    break
                else:
                    prev_para = seq_count

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
        
        FASTA_DIR = [ example_infile_fpath+".split" ]
        if os.path.exists(example_infile_fpath_aligned+".split"):
            FASTA_DIR.append(example_infile_fpath_aligned+".split")

        DIRdict = partition.partition_fasta(
            FASTA_DIR,
            NUMFILE,
            NODESDIR,
            WD,
            WD+"/PARTITION/partition"+str(min(i,MAX_ITERATION-1))+".out",
            "PARTITION.info",
            "UPSTREAM.nwk",
            renamed_subsamplefile_path+".gz",
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
        infile_path,
        infile_path+".aligned",
        iterationfile_path,
        WD+"/tmp.nwk",
        WD+"/SUBSAMPLE.fa.aligned.gz",
        infile_path+".gz.aligned",
        infile_path+".gz.aligned.tree",
        WD+"/seqname_dirpath.txt",
        WD+"/INPUT.part*",
        WD+"/root.fa"
        ]
    
    dirnames = [
        infile_path+".split",
        "EPANG",
        "PARAM",
        "PARTITION",
        "SUBSAMPLE",
        #"TREE"
        ]

    
    subprocess.call(
        "rm -r " + " ".join(filenames+dirnames) + " &> /dev/null",
        shell = True
        )
    
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