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
import manage_edits

import math
import time
import random

def FRACluster(ARGVS, WD, MAX_ITERATION, SUBSAMPLE_SIZE, NODESDIR, THRESHOLD, THREAD_NUM, NUMFILE, QSUBDIR, CODEDIR, 
               ROOTING, MODEL, OPTION,TREEMETHOD, ALIGNED, EPANG, RAXMLSEQ, RAXMLPAR, SOFTWARE,NODE_COUNT,
               INIT_SEQ_COUNT,SEED,ML_or_MP, EXTRACTION_SIZE,careful,FASTA_or_EDIT,
               ALIGNER="unspecified", HMM_PROFILER="unspecified", HMM_ALIGNER="unspecified",
               seq_count_when_aligned=None,
               ):

    print ("FASTA_or_EDIT:",FASTA_or_EDIT)

    # start timer
    start = time.time() 

    ##### hyper parameter #####
    ALIGNMENT_TIMING_PARAMETER = 0.5
    SPLIT_THRESHOLD            = 10000
    mem_req_threshold          = 10**7
    ###########################

    # Enumerate input files
    infile_namelist              = list(sorted(os.listdir(WD)))
    infile_pathlist              = []
    countfile_pathlist           = []
    for infilename in infile_namelist:
        extentions = {FASTA_or_EDIT  , "gz"}
        if infilename.split(".")[-1] in extentions:
            if (infilename!='root.'+FASTA_or_EDIT):
                infile_pathlist.append(WD+"/"+infilename)
        elif infilename.split(".")[-1] == 'count':
            countfile_pathlist.append(WD+"/"+infilename)

    # Record root.fa existed or not
    if (FASTA_or_EDIT == "fa"):
        root_in_separated_file = "root.fa" in set(infile_namelist)
    
    # Create file2Nseq file
    if (len(countfile_pathlist)>0):
        with open(WD + "/file2Nseq.txt", 'w') as ohandle:
            for countfile in countfile_pathlist:
                with open(countfile,'r') as ihandle:
                    ohandle.write(ihandle.read())
                os.remove(countfile)
    # Create file2Nseq dictionary
    if (os.path.exists(WD + "/file2Nseq.txt")):
        print("skip reading files")
        fpath2seqcount = {}
        with open(WD + "/file2Nseq.txt", 'r') as handle:
            for line in handle:
                fpath2seqcount[line.split()[0]] = int(line.split()[1])
    else:
        fpath2seqcount = rename_sequence.count_sequence_fast(infile_pathlist, form = FASTA_or_EDIT)

    ### get input file name ###
    example_infile_fpath         = infile_pathlist[0]
    iterationfile_path           = "unspecified"
    if FASTA_or_EDIT == "fa":
        infile_pathlist_aligned      = [infile_path+".aligned" for infile_path in infile_pathlist]
        example_infile_fpath_aligned = infile_pathlist_aligned[0]
        if (not os.path.isfile(WD+"/root.fa")):
            root_fpath = rename_sequence.outgroup_check_fast(infile_pathlist, "fasta")
        else:
            root_fpath = WD+"/root.fa"
    else:
        example_infile_fpath_aligned = infile_pathlist[0]
        root_fpath = WD+"/root.edit"
    ###########################
    
    ## check input file property ##
    seq_count                 = sum(fpath2seqcount.values())
    is_gzipped                = (example_infile_fpath.split(".")[-1] == "gz")
    if (is_gzipped):
        gzip_extention        = ".gz"
        gzip_command          = "| gzip"
        gunzip_command        = "| gunzip"
    else:
        gzip_extention        = ""
        gzip_command          = ""
        gunzip_command        = ""
    if(INIT_SEQ_COUNT==0): 
        INIT_SEQ_COUNT        = seq_count # only in d0
        seq_count_when_aligned= seq_count
    ################################

    ######## parameter ########
    if(SEED=="random"): random.seed(int(random.randint(0,99999)))
    elif(len(SEED)!=0): random.seed(int(SEED))
    else:print("-r Error: invalid random seed!")
    raxml_thread_num          = 1
    depth                     = max(math.floor(math.log2(seq_count/THRESHOLD))+2,2)
    tree_thread_num           = THREAD_NUM
    ###########################

    '''
    ####### check if aligned #######
    if (ALIGNED=="unaligned"):
        if(os.path.isfile(example_infile_fpath_aligned)):
            if (seq_count < seq_count_when_aligned * ALIGNMENT_TIMING_PARAMETER):
                print(WD,seq_count, seq_count_when_aligned, "alignment needed!!")
                for infile_aligned in infile_pathlist_aligned:
                    os.remove(infile_aligned)
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
        INPUT_FA  = infile_path
    ################################
    '''

    # move to Working Directory
    os.chdir(WD) 
    i=1
    Niteration=1
    
    # call direct tree reconstruction

    if(seq_count<=THRESHOLD):

        concat_infpath = WD+"/INPUT.terminal.fa"

        # Convert .edit(.gz) into .fa
        if (FASTA_or_EDIT == "edit"):
            concat_editpath = WD+"/INPUT.terminal.edit"
            subprocess.call(
                "(cat root.edit "+" ".join(infile_pathlist) + gunzip_command+") > " + concat_editpath,
                shell=True
            )
            edit_list = manage_edits.edit2editlist(concat_editpath)
            manage_edits.edit2fasta(concat_editpath, concat_infpath, edit_list)
            
        elif(FASTA_or_EDIT == "fa"):
            if (root_in_separated_file):
                subprocess.call(
                    "(cat root.fa; cat "+" ".join(infile_pathlist) + gunzip_command+") > "+concat_infpath,
                    shell=True
                )
            else:
                subprocess.call(
                    "(cat "+" ".join(infile_pathlist) + gunzip_command+") > "+concat_infpath,
                    shell=True
                )
        
        if(seq_count<3):
            partition.tiny_tree(concat_infpath,"TERMINAL.nwk")
            print("seq_count < 3!")
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
    
    # call fractal FRACTAL 
    elif(NODE_COUNT>1 and seq_count<=INIT_SEQ_COUNT//NODE_COUNT):

        # To Do: remove concatenation
        if (FASTA_or_EDIT == "fa"):
            concat_infpath = WD+"/INPUT.fa.gz"
            if (root_in_separated_file):
                subprocess.call(
                    "(cat root.fa; cat "+" ".join(infile_pathlist) + gunzip_command+")|gzip > "+concat_infpath,
                    shell=True
                ) 
            else:
                subprocess.call(
                    "(cat "+" ".join(infile_pathlist) + gunzip_command + ")|gzip > " + concat_infpath,
                    shell=True
                )
        elif (FASTA_or_EDIT == "edit"):
            concat_infpath = WD+"/INPUT.edit.gz"
            subprocess.call(
                "(cat "+" ".join(infile_pathlist) + gunzip_command + ")|gzip > " + concat_infpath,
                shell=True
            )

        # Quit distribution after the available computer node saturated 
        os.mkdir("TREE")
        os.chdir(WD+"/TREE")
        if (is_gzipped):
            gzip_option = " -g "
        else:
            gzip_option = ""
        
        open(WD+"/file2Nseq.txt",'w').write(concat_infpath+"\t"+str(seq_count)+"\n")

        FRACTAL_COMMAND =   "FRACTAL"                         + \
                            " -i " + concat_infpath           + \
                            " -k " + str(SUBSAMPLE_SIZE)      + \
                            " -b " + MODEL                    + \
                            " -p " + ML_or_MP                 + \
                            " -P "  + str(careful)            + \
                            " -x " + str(MAX_ITERATION)       + \
                            " -t " + str(THRESHOLD)           + \
                            " -c " + str(THREAD_NUM)          + \
                            " -e "                            + \
                            " -r " + SEED                     + \
                            " -z " + str(EXTRACTION_SIZE)     + \
                            " -n " + WD+"/file2Nseq.txt"       + \
                            gzip_option
        
        if (TREEMETHOD!="unspecified"): 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -m "+TREEMETHOD
        else: 
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -s "+SOFTWARE
        if (ALIGNED=='unaligned'):
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -u "
        if (FASTA_or_EDIT == "edit"):
            FRACTAL_COMMAND = FRACTAL_COMMAND+" -E "
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
        nodenum   = (NODE_COUNT * seq_count) // INIT_SEQ_COUNT - 1

        while i<MAX_ITERATION:

            Niteration = i
            
            os.chdir(WD)

            ### set file names###
            iterationfile_path         = WD + "/ITERATION." + FASTA_or_EDIT                   
            subsamplefile_path         = WD + "/SUBSAMPLE/SUBSAMPLE." + FASTA_or_EDIT             
            renamed_subsamplefile_path = WD + "/SUBSAMPLE/RENAMED_"+str(i)+"." + FASTA_or_EDIT    
            #####################
        
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
            if (FASTA_or_EDIT == "fa"):
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
                        if fpath2seqcount[file_path] > Nseq_per_file:
                            if (FASTA_or_EDIT == 'fa'):
                                subprocess.call("seqkit split2 -s "+str(Nseq_per_file)+" "+file_path+" &> /dev/null; rm "+file_path, shell=True)
                            elif (FASTA_or_EDIT == 'edit'):
                                os.mkdir(file_path+".split")
                                subprocess.call(
                                    "cd " +file_path + ".split; "+ 
                                    "cat "+file_path + gunzip_command+"|split -l "+str(Nseq_per_file)+" &> /dev/null;"+
                                    "rm " +file_path +
                                    "for file in $(ls); do cat $file "+gzip_command+" >"+file_path.split(".")[-1]+".${file}"+gzip_extention+"; done", 
                                    shell=True)
                        else:
                            os.mkdir(file_path+".split")
                            shutil.move(file_path, file_path+".split/")
                        splitted_fnamelist = sorted(os.listdir(file_path+".split"))
                        for k, splitted_fname in enumerate(splitted_fnamelist):
                            if (k < len(splitted_fnamelist)-1):
                                fpath2seqcount[splitted_dirpath + "/" + splitted_fname] = Nseq_per_file
                            else:
                                if k > 0:
                                    fpath2seqcount[splitted_dirpath + "/" + splitted_fname] = fpath2seqcount[file_path] % Nseq_per_file
                                else:
                                    fpath2seqcount[splitted_dirpath + "/" + splitted_fname] = fpath2seqcount[file_path]
                        fpath2seqcount.pop(file_path)
                    for file_path in file_pathlist_to_be_splitted[1:]:
                        subprocess.call(
                            "for file in $(ls "+ file_path+".split/); do mv "+ file_path+".split/$file " + splitted_dirpath + "; done; "
                            "rm -r " + file_path+".split",
                            shell=True
                        )
            else:
                if (not os.path.exists(example_infile_fpath+".split")):
                    os.mkdir(example_infile_fpath+".split")
                    for infile_path in infile_pathlist:
                        fpath2seqcount[example_infile_fpath+".split/"+infile_path.split("/")[-1]] = fpath2seqcount[infile_path]
                        fpath2seqcount.pop(infile_path)
                        shutil.move(infile_path, example_infile_fpath+".split")
                    splitted_dirpath = example_infile_fpath+".split"
            
            # rename
            
            for j, splitted_fname in enumerate(os.listdir(splitted_dirpath)):
                splitted_fpath   = splitted_dirpath+"/"+splitted_fname
                renamed_filepath = splitted_dirpath+"/INPUT.part"+str(j)+"."+WD.split("/")[-1]+"."+str(i)+"."+FASTA_or_EDIT+gzip_extention
                if renamed_filepath != splitted_fpath:
                    shutil.move(splitted_fpath, renamed_filepath)
                    fpath2seqcount[renamed_filepath] = fpath2seqcount[splitted_fpath]
                    fpath2seqcount.pop(splitted_fpath)
            

            #################
            #random sampling#
            #################

            if (os.path.exists(iterationfile_path)):
                sampled_seq_name_list = \
                    rename_sequence.random_sampling( # fasta or edit
                        in_dirpath = None,
                        out_fname  = subsamplefile_path,
                        subsample_size = SUBSAMPLE_SIZE,
                        fpath2seqcount = None,
                        root_fpath     = root_fpath,
                        total_seqcount = None, 
                        file_format    = FASTA_or_EDIT,
                        in_fpath       = iterationfile_path
                        )
            else:
                # subsampling
                sampled_seq_name_list = \
                    rename_sequence.random_sampling( # fasta or edit
                        in_dirpath = splitted_dirpath,
                        out_fname  = subsamplefile_path,
                        subsample_size = SUBSAMPLE_SIZE,
                        fpath2seqcount = fpath2seqcount,
                        root_fpath     = root_fpath,
                        total_seqcount = seq_count, 
                        file_format    = FASTA_or_EDIT,
                        in_fpath       = None
                        )
            
            #################
            #rename sequence#
            #################
            seqname2renamedname = rename_sequence.rename_sequence(
                subsamplefile_path, 
                renamed_subsamplefile_path,
                file_format = FASTA_or_EDIT,
            )

            if (FASTA_or_EDIT == 'edit'):
                edit_list = manage_edits.edit2editlist(renamed_subsamplefile_path)
                manage_edits.edit2fasta(renamed_subsamplefile_path, renamed_subsamplefile_path+".fa", edit_list)
                renamed_subsamplefile_path = renamed_subsamplefile_path+".fa"
            
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
                
                extracted_seq_name_list   = ["root"]+random.sample(sampled_seq_name_list, EXTRACTION_SIZE-1)
                extracted_seq_rename_list = [ seqname2renamedname[name] for name in extracted_seq_name_list ]
                extraction.tree_extraction(
                    renamed_subsamplefile_path+".aligned.tree",
                    set(extracted_seq_rename_list)                         ,
                    renamed_subsamplefile_path+".aligned.extracted.tree"
                    )
                extraction.fasta_extraction(
                    renamed_subsamplefile_path          ,
                    set(extracted_seq_rename_list)                      ,
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted",
                    )
                extraction.fasta_extraction(
                    renamed_subsamplefile_path+".aligned"          ,
                    set(extracted_seq_rename_list)                      ,
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.aligned",
                    )
                os.remove(renamed_subsamplefile_path+".aligned.tree")
                os.remove(renamed_subsamplefile_path+"")
                os.remove(renamed_subsamplefile_path+".aligned")
                
                os.rename(
                    renamed_subsamplefile_path+".aligned.extracted.tree",
                    renamed_subsamplefile_path+".aligned.tree"
                    )
                shutil.copyfile(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted",
                    renamed_subsamplefile_path
                )
                os.rename(
                    WD + "/SUBSAMPLE/RENAMED_"+str(i)+".fa.extracted.aligned",
                    renamed_subsamplefile_path+".aligned"
                )

            os.chdir(WD+"/PARAM")
            
            if(raxml_thread_num>1):
                subprocess.call(
                    RAXMLPAR                                                        +
                    " -T "   + str(raxml_thread_num)                                +
                    " -f e"                                                         +
                    " -s "   + renamed_subsamplefile_path+".aligned"                +
                    " -t "   + renamed_subsamplefile_path+".aligned.tree"           +
                    " -n "   + "PARAM_"+str(i)                                      +
                    " -m "   + MODEL                                                +
                    " 2> /dev/null"                                                ,
                    shell=True
                )
            else:
                subprocess.call(
                    RAXMLSEQ                                                        +
                    " -f e"                                                         +
                    " -s "   + renamed_subsamplefile_path+".aligned"                + 
                    " -t "   + renamed_subsamplefile_path+".aligned.tree"           +
                    " -n "   + "PARAM_"+str(i)                                      +
                    # " -n "   + "PARAM_"+str(i-i%2)                                  + # for test
                    " -m "   + MODEL                                                +
                    " 2> /dev/null"                                                 ,
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
                if(os.path.exists(example_infile_fpath_aligned+".split")):
                    QUERY_FA_DIR          = example_infile_fpath_aligned+".split"
                    ALIGNED_FOR_PLACEMENT = "aligned"
                else:
                    QUERY_FA_DIR          = example_infile_fpath+".split"
                    ALIGNED_FOR_PLACEMENT = ALIGNED
                
                if (ALIGNED=="aligned"):
                    refseq = renamed_subsamplefile_path
                else:
                    refseq = renamed_subsamplefile_path+".aligned"

                if FASTA_or_EDIT == 'fa':
                    edit_list = None

                # conduct placement
                placemnt_exit_code = placement.distributed_placement(
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
                    file_format = FASTA_or_EDIT                         ,
                    edit_list = edit_list                               ,
                    alignment_outdir = example_infile_fpath_aligned+".split",
                    careful=careful                                     ,
                    hmm_aligner=HMM_ALIGNER                             ,
                    hmm_profiler=HMM_PROFILER
                )

                if (placemnt_exit_code == 1):

                    print("Placement failed...")

                    # if i == 0, start from random sampling again, else use the result of previous i
                    if(i > 1 and os.path.isfile(WD+"/PARTITION/partition"+str(i-1)+".out") ): 
                        i -= 1
                        para  = prev_para
                        break
                    else:
                        if(os.path.isfile(iterationfile_path)):
                            os.remove(iterationfile_path)
                        i += 1
                
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
                if(para>prev_para or Nseq_in_largest_subclade == seq_count):
                    
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

                if(para!=0): # if problematic sequences remained
                    
                    if (FASTA_or_EDIT == 'fa'):
                        # select subsample sequence file
                        if os.path.isfile(subsamplefile_path+".aligned"):
                            ALIGNED_SUBSAMPLE = subsamplefile_path+".aligned"
                        else:
                            ALIGNED_SUBSAMPLE = subsamplefile_path
                            
                        subprocess.call(
                            "cat "+ALIGNED_SUBSAMPLE+" "+WD+"/EPANG/problematic.fa > "+iterationfile_path,
                            shell=True
                        )
                        
                        # select all sequence file
                        if os.path.exists(example_infile_fpath_aligned+".split"):
                            ALIGNED_ALL_DIR = example_infile_fpath_aligned+".split"
                            ALIGNED         = "aligned"
                        else:
                            ALIGNED_ALL_DIR = example_infile_fpath        +".split"

                    elif (FASTA_or_EDIT == 'edit'):
                        shutil.copyfile(
                            WD+"/SUBSAMPLE/RENAMED_"+str(i)+".edit",
                            WD+"/ITERATION.edit"
                        )

                        partition.add_paraphyletic_edit(
                            WD+"/PARTITION/partition"+str(i)+".out",
                            "ITERATION.edit"                       ,
                            splitted_dirpath                       ,
                            SUBSAMPLE_SIZE                         ,
                            para                                   ,
                            )

                    i+=1
                    prev_para=para
                elif (not resampling_needed):
                    break
                else:
                    prev_para = seq_count

        ###########################
        #partition .fasta or .edit#
        ###########################
        os.chdir(WD)
        try:
            print(
                "<Sequence count> Input:" + str(seq_count)                + "\n" +
                "  Largest subclade    :" + str(Nseq_in_largest_subclade) + "\n" +
                "  Problematic         :" + str(para)
            )
        except:
            None
        
        if (Nseq_in_largest_subclade == seq_count ): # if all sequences were classified into one subclade, FRACTAL gives up for inference of this clade
            print ("Error: FRACluster.py could not divide sequences into multiple subclades")
            return
        
        FASTA_DIR = [ example_infile_fpath+".split" ]
        if os.path.exists(example_infile_fpath_aligned+".split") and FASTA_or_EDIT=='fa':
            FASTA_DIR.append(example_infile_fpath_aligned+".split")

        DIRdict = partition.partition_fasta(
            FASTA_DIR,
            NUMFILE,
            NODESDIR,
            WD,
            WD+"/PARTITION/partition"+str(min(i,MAX_ITERATION-1))+".out",
            "UPSTREAM.nwk",
            ROOTING,
            nodenum = nodenum,
            codedir = CODEDIR,
            file_format = FASTA_or_EDIT
            )
        
        dirpath2Nseq_filepath = WD + "/Nseq_dirpath.txt"

        partition.qsub_prep(
            ARGVS,
            WD, 
            DIRdict, 
            INIT_SEQ_COUNT,
            seq_count_when_aligned,
            dirpath2Nseq_filepath,
            mem_req_threshold
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
        example_infile_fpath,
        example_infile_fpath+".aligned",
        iterationfile_path,
        WD+"/tmp.nwk",
        WD+"/SUBSAMPLE.fa.aligned",
        WD+"/SUBSAMPLE.fa.aligned.gz",
        example_infile_fpath+".gz.aligned",
        example_infile_fpath+".gz.aligned.tree",
        WD+"/seqname_dirpath.txt",
        WD+"/INPUT.part*",
        WD+"/INPUT.terminal*",
        WD+"/root.fa",
        WD+"/INPUT.fa",
        WD+"/INPUT.fa.gz",
        WD+"/INPUT.fa*aligned*"
        ]
    
    dirnames = [
        example_infile_fpath+".split",
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
            str(Niteration)   + "," +
            "subsamplings"    + "\n"
        )

if __name__ == "__main__":
    argvs = sys.argv

    try:
        sys.setrecursionlimit(1000000000)
    except:
        None

    print(argvs)

    if (len(argvs)==27):
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
            FASTA_or_EDIT     = argvs[25],
            ALIGNER           = "unspecified", 
            HMM_PROFILER      = "unspecified", 
            HMM_ALIGNER       = "unspecified",
            seq_count_when_aligned=None
            )
    elif ((len(argvs)==27+3)):
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
            FASTA_or_EDIT     = argvs[25],
            ALIGNER           = argvs[26],
            HMM_PROFILER      = argvs[27],
            HMM_ALIGNER       = argvs[28],
            seq_count_when_aligned=int(argvs[29])
        )
    else:
        print("Error: Number of arguments: "+str(len(argvs))+" for FRACluster.py is wrong!")