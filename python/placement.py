from Bio import SeqIO
import jplace_parse
import divide_ref_and_query
import manage_edits
import subprocess
import os
import shutil
import gzip

# decompose fasta file (in_file) into x fasta files.
def decompose_fasta(in_file, x,seq_count):
    n=seq_count
    k=n//x # each decomposed file has k sequences
    ohandle=[]
    for i in range(x):
        ohandle.append(
            gzip.open(in_file+"."+str(i)+".gz", 'wt')
        ) 
    with gzip.open(in_file,'rt') as ihandle:
        allseq_itr = SeqIO.parse(ihandle, "fa")
        l=0 # inclement constantly
        if(True):
            for record in allseq_itr:
                i = min(l//k, x-1)
                SeqIO.write(record,ohandle[i],'fa')
                l+=1
    for i in range(x):
        ohandle[i].close()

def decompose_edit(in_file, x,seq_count):
    n=seq_count
    k=n//x # each decomposed file has k sequences
    ohandle=[]
    for i in range(x):
        ohandle.append(
            gzip.open(in_file+"."+str(i)+".gz", 'wt')
        ) 
    with gzip.open(in_file,'rt') as ihandle:
        l=0 # inclement constantly
        if(True):
            for line in ihandle:
                i = min(l//k, x-1)
                ohandle[i].write(line)
                l+=1
    for i in range(x):
        ohandle[i].close()

def decompose_edit2(in_file, seq_count, n_per_file = 10000):
    n=seq_count
    k=10000
    x=n//k
    ohandle=[]
    for i in range(x):
        ohandle.append(
            open(in_file+"."+str(i)+".gz", 'w')
        ) 
    with open(in_file,'r') as ihandle:
        l=0 # inclement constantly
        if(True):
            for line in ihandle:
                i = min(l//k, x-1)
                ohandle[i].write(line)
                l+=1
    for i in range(x):
        ohandle[i].close()

def distributed_placement(  WD, EPANG, refseq, reftree, model, 
                            query_dir, outdir, threadnum, nodenum, 
                            codedir, seq_count, ML_or_MP, RAXMLSEQ, 
                            ALIGNED, seed, careful=1, hmm_aligner="", hmm_profiler="",
                            file_format = "fa", edit_list = None, alignment_outdir = None):
    
    splitted_queryfile_list = os.listdir(query_dir)

    # Profile HMM
    if(ALIGNED=="unaligned"): # for unaligned sequences
        # Build HMM profile
        subprocess.call(
            hmm_profiler + " " + 
            refseq+".hmm "     +
            refseq             ,
            shell=True
        ) 
        os.mkdir(alignment_outdir)

    # sequential mode
    if(nodenum<=1):

        for filename in splitted_queryfile_list:

            query =  query_dir + "/" + filename
            os.mkdir(outdir    + "/" + filename)
            os.chdir(outdir    + "/" + filename)
            if query.split(".")[-1] == "gz":
                is_gzipped = True
                extention  = ".gz"
                gzipcommand= "|gzip"
                gunzipcommand = "|gunzip"
            else:
                is_gzipped = False
                extention  = ""
                gzipcommand= ""
                gunzipcommand = ""

            if (file_format == "edit"):
                manage_edits.edit2fasta(query, outdir + "/" + filename + "/" + filename + ".fa.gz", edit_list) # TO DO
                query = outdir + "/" + filename + "/" + filename + ".fa.gz"

            if(ML_or_MP=="ML"): 
                if(ALIGNED=="unaligned"):
                    # One-by-one HMM alignment
                    print(
                        hmm_aligner+" "        +
                        "--outformat afa "     +
                        "--mapali "+refseq+" " +
                        "--trim "              + # for trimming insersions
                        refseq+".hmm "         +
                        query +" "             +
                        "| sed 's/\./-/g'> "   +
                        outdir+"/"+filename+"/ref_query.fa"
                    )
                    subprocess.call(
                        hmm_aligner+" "        +
                        "--outformat afa "     +
                        "--mapali "+refseq+" " +
                        "--trim "              + # for trimming insersions
                        refseq+".hmm "         +
                        query +" "             +
                        "| sed 's/\./-/g'> "   +
                        outdir+"/"+filename+"/ref_query.fa"
                        ,shell=True
                    )
                    divide_ref_and_query.\
                        divide_fasta_into_ref_and_query(
                            outdir+"/"+filename+"/ref_query.fa", 
                            refseq
                        )
                    subprocess.call(
                        EPANG                                            +
                        " --redo"                                        +
                        " -s "+outdir+"/"+filename+"/ref_query.fa.ref"   +
                        " -t "+reftree                                   +
                        " --model "+model                                +
                        " -q "+outdir+"/"+filename+"/ref_query.fa.query" +
                        " -w "+outdir+"/"+filename                       +
                        " -T "+str(threadnum)                            ,
                        shell=True
                    )
                elif(ALIGNED=="aligned"): # for aligned sequences
                    subprocess.call(
                        EPANG                  +
                        " --redo"              +
                        " -s "+refseq          +
                        " -t "+reftree         +
                        " --model "+model      +
                        " -q "+query           +
                        " -w "+outdir+"/"+filename         +
                        " -T "+str(threadnum)  ,
                        shell=True
                    )
                jplace_parse.parse_jplace(
                    outdir+"/"+filename+"/epa_result.jplace",
                    "epa-ng",
                    query,
                    seed,
                    careful=careful
                    )
                
            if(ML_or_MP=="MP"): 
                if(ALIGNED=="unaligned"): # for unaligned sequences
                    # One-by-one HMM alignment
                    subprocess.call(
                        hmm_aligner                +
                        " --outformat afa"         +
                        " --mapali "+refseq+" "    +
                        refseq+".hmm "             +
                        query                      +
                        " | sed 's/\./N/g'> "      +
                        outdir+"/"+filename+"/ref_query.fa"  ,
                        shell=True
                    )  
                elif(ALIGNED=="aligned"): # for aligned sequences
                    subprocess.call(
                        "(cat "+refseq+"; cat "      +
                        query + gunzipcommand +")"   +
                        " > "+outdir+"/"+filename+"/ref_query.fa" ,
                        shell=True
                    )
                os.chdir(outdir+"/"+filename)
                subprocess.call(
                    RAXMLSEQ                      +
                    " -n epa_result"              +
                    " -f y -m GTRCAT"             +
                    " -s "+outdir+"/"+filename+"/ref_query.fa" +
                    " -t "+reftree                ,
                    shell=True
                )
                jplace_parse.parse_jplace(
                    outdir+"/"+filename+"/RAxML_portableTree.epa_result.jplace",
                    "epa_MP",
                    outdir+"/"+filename+"/ref_query.fa",
                    seed,
                    careful=careful
                )
            # If HMM alignments were conducted
            if(ALIGNED=="unaligned"):
                subprocess.call(
                    "cat " + outdir+"/"+filename+"/ref_query.fa.query" + gzipcommand + ">" + 
                    alignment_outdir+"/"+filename+".aligned"+extention,
                    shell = True
                    )
                shutil.move(
                    outdir+"/"+filename+"/ref_query.fa.ref", 
                    WD+"/SUBSAMPLE.fa.aligned"
                    )
        
        subprocess.call(
            "for dir in $(ls "+outdir+"); do "+
            "   mv "  + outdir + "/$\{dir\}/placement_tree.out "+outdir+"/placement_tree.out;"+
            "   cat " + outdir + "/$\{dir\}/edge_to_seqname.out >> " + outdir+"/edge_to_seqname_all.out;"        +
            "   cat " + outdir + "/$\{dir\}/problematic.fa >> " + outdir+"/problematic.fa;"+
            "done",
            shell=True
            )

    else: # in distributed computing mode
        dname=WD.split("/").pop()
        
        if ( file_format == "fa" ):
            #moved=outdir+"/query.fa.gz"
            #shutil.move(query, moved)

            if(ALIGNED=="unaligned"):
                # Build HMM profile
                subprocess.call(
                    hmm_profiler+" "+
                    refseq+".hmm "+
                    refseq,
                    shell=True
                    )
        elif ( file_format == "edit" ):
            moved=outdir+"/query.edit.gz"
            shutil.move(query, moved) # TODO
            decompose_edit(
                moved,
                nodenum, 
                seq_count
            )
            node2filelist = [ ["query.edit.gz."+str(i)+".gz"] for i in range(nodenum) ]
            splitted_files_dir = outdir + "/"
            with open(outdir+"/editlist.txt", 'w') as handle:
                for edit in edit_list:
                    handle.write(edit + "\n")

        Nfiles_total = len(splitted_queryfile_list)
        Nfiles_per_node = len(splitted_queryfile_list) // nodenum # Only the last node may treat more number of files
        node2filelist = []
        for i in range(nodenum):
            node2filelist.append([])
        for j, file_name in enumerate(splitted_queryfile_list):
            node2filelist[j%nodenum].append(file_name)

        #distribution start
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
                Popen(  'echo $LD_LIBRARY_PATH', 
                        stdout=subprocess.PIPE,
                        shell=True
                ).communicate()[0]
            ).decode('utf-8')
        LD_LIBRARY_PATH = (LD_LIBRARY_PATH.split('\n'))[0]
        for i in range(nodenum):
            os.mkdir(outdir+"/EPANG"+str(i))
            with open(WD+"/../../prep_dir/qsub_"+dname+"."+str(i)+".placement.sh", 'w') as handle:
                for filename in node2filelist[i]:

                    if filename.split(".")[-1] == "gz":
                        is_gzipped = True
                        extention  = ".gz"
                        gzipcommand= "|gzip"
                        gunzipcommand = "|gunzip"
                    else:
                        is_gzipped = False
                        extention  = ""
                        gzipcommand= ""
                        gunzipcommand = ""


                    os.mkdir(outdir+"/EPANG"+str(i)+"/"+filename)
                    queryfile = query_dir + "/" + filename

                    handle.write("#!/bin/bash\n")
                    handle.write("#$ -S /bin/bash\n")
                    handle.write("PATH={}\n".format(PATH))
                    handle.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))

                    if (file_format == "edit"):
                        handle.write(
                            "python3 "                              +
                            codedir+"/python/manage_edits.py "      +
                            queryfile + " "                         +
                            outdir +"/editlist.txt\n"
                            )
                        queryfile = outdir+"/EPANG"+str(i)+"/"+filename+"/query.edit.fa.gz" + "." + str(i) + ".gz"
                        handle.write(
                            "mv "                                        + \
                            outdir+"/query.edit.gz."+str(i)+".gz.fa.gz " + \
                            queryfile + "\n"
                            )

                    if(ML_or_MP=="ML"): 
                        if(ALIGNED=="unaligned"): # for unaligned sequences
                            # Conduct HMM alignment
                            handle.write(
                                hmm_aligner                                +
                                " --outformat afa"                         +
                                " --mapali " + refseq + " "                +
                                refseq+".hmm "                             +
                                queryfile                                  +
                                " | sed 's/\./-/g' "                       +
                                " > "                                      +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa\n"
                                )
                            handle.write(
                                "python3 "                                 +
                                codedir+"/python/divide_ref_and_query.py " +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa "    + 
                                refseq + "\n"
                                )
                            handle.write(
                                "trimal "                                                  +
                                " -in " + outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa "       +
                                "-selectcols { "                                           +
                                "   `trimal -sgc "                                         +
                                "    -in " +outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.ref " +
                                "    |awk ' { if( $2==100 ){ print $1 }}'"                 +
                                "    |tr \"\\n\" \",\" | sed -e \"s/,\$//\" ` "            +
                                "    } "                                                   +
                                "> "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols\n"
                                )
                            handle.write(
                                "python3 "                                          +
                                codedir +"/python/divide_ref_and_query.py "         +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols "+
                                refseq  + "\n"
                                )
                            
                            handle.write(
                                EPANG                                               + 
                                " --redo"                                           +
                                " -s "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.ref"   +
                                " -t "+reftree                                      +
                                " --model "+model                                   +
                                " -q "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.query" +
                                " -w "+outdir+"/EPANG"+str(i)+"/"+filename                       +
                                " -T "+str(threadnum)+"\n"
                                )
                            handle.write(
                                "cat " + outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols.query" + gzipcommand +
                                ">" + 
                                alignment_outdir+"/"+filename+".aligned"+extention+"\n"
                                )
                        elif(ALIGNED=="aligned"): # for aligned sequences
                            handle.write(
                                EPANG                         +
                                " --redo"                     +
                                " -s "+refseq                 +
                                " -t "+reftree                +
                                " --model "+model             +
                                " -q "+queryfile              +
                                " -w "+outdir+"/EPANG"+str(i)+"/"+filename +
                                " -T "+str(threadnum)+"\n"
                                )
                        handle.write(
                            "cd "+outdir+"/EPANG"+str(i)+"/"+filename+"\n"
                            )
                        handle.write(
                            "python3 "                                   +
                            codedir+"/python/jplace_parse.py "           +
                            outdir+"/EPANG"+str(i)+"/"+filename+"/epa_result.jplace " +
                            "epa-ng "                                    +
                            queryfile + " "                              + 
                            seed                                         +
                            "\n"
                            )
                    elif(ML_or_MP=="MP"):
                        handle.write(
                            "cd "+outdir+"/EPANG"+str(i)+"/"+filename+"\n"
                        )
                        if(ALIGNED=="unaligned"): # for unaligned sequences
                            # Conduct HMM alignment
                            handle.write(
                                hmm_aligner                                 +
                                " --outformat afa"                          +
                                " --mapali "+refseq+" "                     +
                                refseq+".hmm "+queryfile                    +
                                " | sed 's/\./-/g' "                        +
                                " > "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa\n"
                            )
                            handle.write(
                                "cat "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa " +
                                "| sed 's/\./N/g' > " + outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.gap2N.fa\n"
                            )
                            handle.write(
                                "python3 "                                 +
                                codedir+"/python/divide_ref_and_query.py " +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa "    + 
                                refseq + "\n"
                                )
                            handle.write(
                                "trimal "                                                  +
                                " -in " + outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa "       +
                                "-selectcols { "                                           +
                                "   `trimal -sgc "                                         +
                                "    -in " +outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.ref " +
                                "    |awk ' { if( $2==100 ){ print $1 }}'"                 +
                                "    |tr \"\\n\" \",\" | sed -e \"s/,\$//\" ` "            +
                                "    } "                                                   +
                                "> "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols\n"
                                )
                            handle.write(
                                "python3 "                                          +
                                codedir +"/python/divide_ref_and_query.py "         +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols "+
                                refseq  + "\n"
                                )
                        elif(ALIGNED=="aligned"): # for aligned sequences
                            handle.write(
                                "(cat "+refseq+"; cat "+
                                queryfile        +
                                gunzipcommand    +
                                ")| sed 's/\./N/g'> "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.gap2N.fa\n"
                            )
                        handle.write(
                            RAXMLSEQ                                      +
                            " -n epa_result -f y -m GTRCAT"               +
                            " -s "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.gap2N.fa" +
                            " -t "+reftree+"\n"
                            ) 
                        handle.write(
                            "python3 "                                                      +
                            codedir+"/python/jplace_parse.py "                              +
                            outdir+"/EPANG"+str(i)+"/"+filename+"/RAxML_portableTree.epa_result.jplace " +
                            "epa_MP "                                                       +
                            outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.gap2N.fa "      +
                            seed + "\n"
                            )
                        files_to_be_removed = [
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"RAxML_equallyParsimoniousPlacements.epa_result",
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"RAxML_originalLabelledTree.epa_result"         ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"RAxML_portableTree.epa_result.jplace"          ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"RAxML_labelledTree.epa_result"                 ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"RAxML_labelledTree.epa_result"                 ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"ref_query*fa"                                  ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"*reduced"                                      ,
                            outdir+"/EPANG"+str(i)+"/"+filename+"/"+"*.selectcols"                                  ,
                        ]
                        handle.write(
                            "rm "+ " ".join(files_to_be_removed) + " &> /dev/null\n"
                            ) 
                        if(ALIGNED!="aligned"):
                            handle.write(
                                "cat "                                                               +
                                outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.selectcols.query" +
                                gzipcommand + " > "                                                  +
                                alignment_outdir+"/"+filename+".aligned" + extention + "\n"          ,
                            )
                handle.write(
                    "echo \"finished\" > "      +
                    outdir+"/epang"+str(i)+".o\n"
                    )
                # end of a distributed task
            shutil.move(WD+"/../../prep_dir/qsub_"+dname+"."+str(i)+".placement.sh", WD+"/../../qsub_dir/qsub_"+dname+"."+str(i)+".placement.sh")
        # check if all placement tasks ended
        flag = 0
        while(flag==0):
            i=0
            while(i<nodenum):
                if(not os.path.exists(outdir+"/epang"+str(i)+".o")):
                    break
                i+=1
            if i == nodenum:
                flag=1
        
        # remove unnecessary files
        for i in range(nodenum):
            os.remove(outdir+"/epang"+str(i)+".o")
        
        if (file_format == "fa"):
            None
        elif (file_format == "edit"):
            shutil.move(outdir+"/query.edit.gz",query)

        # (If HMM alignments were conducted) concat aligned sequences
        if(ALIGNED=="unaligned"):
            subprocess.call(
                "cat "+
                "$(ls "+outdir+"/EPANG0/*/ref_query.fa.selectcols.ref | head -n1)"+
                "> "+ 
                WD+"/SUBSAMPLE.fa.aligned",
                shell=True
                )

        # merge results
        subprocess.call(
            "mv "  + "$(ls "+outdir+"/EPANG0/*/placement_tree.out | head -n1) "+outdir+"/placement_tree.out;",
            shell=True
            )
        for i in range(nodenum):
            subprocess.call(
                "cat " + outdir + "/EPANG"+str(i)+"/*/edge_to_seqname.out >> " + outdir+"/edge_to_seqname_all.out; "       +
                "cat " + outdir + "/EPANG"+str(i)+"/*/problematic."+file_format+" >> " + outdir+"/problematic."+file_format,
                shell=True
                )

def my_paste(outdir, nodenum, outfilename):
    handles=[]
    for i in range(nodenum):
        if os.path.exists(outdir+"/EPANG"+str(i)+"/edge_to_seqname.out"):
            handles.append(open(outdir+"/EPANG"+str(i)+"/edge_to_seqname.out"))
    with open(outfilename,'w') as outhandle:
        line=handles[0].readline()
        while line:
            mergedline=line.split("\n")[0]
            for i in range(1,len(handles)):
                line=handles[i].readline().split("\n")[0]
                mergedline+=line
            outhandle.write(mergedline+"\n")
            line=handles[0].readline()
    for i in range(len(handles)):
        handles[i].close()