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
        allseq_itr = SeqIO.parse(ihandle, "fasta")
        l=0 # inclement constantly
        if(True):
            for record in allseq_itr:
                i = min(l//k, x-1)
                SeqIO.write(record,ohandle[i],'fasta')
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
                            query, outdir, threadnum, nodenum, 
                            codedir, seq_count, ML_or_MP, RAXMLSEQ, 
                            ALIGNED, seed, careful=1, hmm_aligner="", hmm_profiler="",
                            file_format = "fasta", edit_list = None):
    if(nodenum<=1):

        if (file_format == "edit"):
            manage_edits.edit2fasta(query, edit_list)
            query = query + ".fa.gz"

        if(ML_or_MP=="ML"): 
            if(ALIGNED=="unaligned"): # for unaligned sequences
                # Build HMM profile
                subprocess.call(
                    hmm_profiler + " " + 
                    refseq+".hmm "     +
                    refseq             ,
                    shell=True
                ) 
                # One-by-one HMM alignment
                subprocess.call(
                    hmm_aligner+" "        +
                    "--outformat afa "     +
                    "--mapali "+refseq+" " +
                    "--trim "              + # for trimming insersions
                    refseq+".hmm "         +
                    query +" "             +
                    "| sed 's/\./-/g'> "   +
                    outdir+"/ref_query.fa"
                    ,shell=True
                )
                divide_ref_and_query.\
                    divide_fasta_into_ref_and_query(
                        outdir+"/ref_query.fa", 
                        refseq
                    )
                subprocess.call(
                    EPANG                               +
                    " --redo"                           +
                    " -s "+outdir+"/ref_query.fa.ref"   +
                    " -t "+reftree                      +
                    " --model "+model                   +
                    " -q "+outdir+"/ref_query.fa.query" +
                    " -w "+outdir                       +
                    " -T "+str(threadnum)               ,
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
                    " -w "+outdir          +
                    " -T "+str(threadnum)  ,
                    shell=True
                )
            os.chdir(outdir)
            jplace_parse.parse_jplace(
                outdir+"/epa_result.jplace",
                "epa-ng",
                seed,
                careful=careful
                )
            
        if(ML_or_MP=="MP"): 
            if(ALIGNED=="unaligned"): # for unaligned sequences
                # Build HMM profile
                subprocess.call(
                    hmm_profiler+" "  +
                    refseq+".hmm "    +
                    refseq            ,
                    shell=True
                )
                # One-by-one HMM alignment
                subprocess.call(
                    hmm_aligner                +
                    " --outformat afa"         +
                    " --mapali "+refseq+" "    +
                    refseq+".hmm "             +
                    query                      +
                    " | sed 's/\./N/g'> "      +
                    outdir+"/ref_query.fa"  ,
                    shell=True
                )  
            elif(ALIGNED=="aligned"): # for aligned sequences
                subprocess.call(
                    "cat "+refseq+" "            +
                    query                        +
                    " | gunzip > "+outdir+"/ref_query.fa" ,
                    shell=True
                )
            os.chdir(outdir)
            subprocess.call(
                RAXMLSEQ                      +
                " -n epa_result"              +
                " -f y -m GTRCAT"             +
                " -s "+outdir+"/ref_query.fa" +
                " -t "+reftree                ,
                shell=True
            )
            jplace_parse.parse_jplace(
                outdir+"/RAxML_portableTree.epa_result.jplace",
                "epa_MP",
                seed,
                careful=careful
            )
        os.rename(
            outdir+"/edge_to_seqname.out",
            outdir+"/edge_to_seqname_all.out"
            )
        # If HMM alignments were conducted
        if(ALIGNED=="unaligned"):
            shutil.move(
                outdir+"/ref_query.fa.query", 
                WD+"/INPUT.fa.aligned"
                )
            shutil.move(
                outdir+"/ref_query.fa.ref", 
                WD+"/SUBSAMPLE.fa.aligned"
                )

    else: # in distributed computing mode
        dname=WD.split("/").pop()
        
        if ( file_format == "fasta" ):
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
            #decompose_fasta(
            #    moved,
            #    nodenum, 
            #    seq_count
            #    )
        elif ( file_format == "edit" ):
            moved=outdir+"/query.edit.gz"
            shutil.move(query, moved)
            decompose_edit(
                moved,
                nodenum, 
                seq_count
            )
            with open(outdir+"/editlist.txt", 'w') as handle:
                for edit in edit_list:
                    handle.write(edit + "\n")
        
        if (os.path.exists(WD + "/INPUT.fa.aligned.gz.split")):
            splitted_fasta_dir  = WD + "/INPUT.fa.aligned.gz.split/"
            splitted_fasta_list = os.listdir(WD + "/INPUT.fa.aligned.gz.split")
            splitted = True
        elif (os.path.exists(WD + "/INPUT.fa.gz.split")):
            splitted_fasta_dir  = WD + "/INPUT.fa.gz.split/"
            splitted_fasta_list = os.listdir(WD + "/INPUT.fa.gz.split")
            splitted = True
        if (splitted):
            Nfiles_total = len(splitted_fasta_list)
            Nfiles_per_node = len(splitted_fasta_list) // nodenum + 1 # Only the last node may treat more number of files
            node2filelist = []
            for i in range(nodenum):
                if   (i <  nodenum - 1):
                    file_list = splitted_fasta_list[Nfiles_per_node * i:Nfiles_per_node * (i+1)]
                elif (i == nodenum - 1):
                    file_list = splitted_fasta_list[Nfiles_per_node * i:]
                node2filelist.append(file_list)

        #distribution start
        for i in range(nodenum):
            os.mkdir(outdir+"/EPANG"+str(i))
            with open(WD+"/../../qsub_dir/qsub_"+dname+"."+str(i)+".sh", 'w') as handle:
                for filename in node2filelist[i]:
                    os.mkdir(outdir+"/EPANG"+str(i)+"/"+filename)
                    queryfile = splitted_fasta_dir + filename

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
                    handle.write("#!/bin/bash\n")
                    handle.write("#$ -S /bin/bash\n")
                    handle.write("PATH={}\n".format(PATH))
                    handle.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))

                    if (file_format == "edit"):
                        handle.write(
                            "python3 "                              +
                            codedir+"/python/manage_edits.py "      +
                            outdir +"/query.edit.gz."+str(i)+".gz " +
                            outdir +"/editlist.txt\n"
                            )
                        moved = outdir+"/query.edit.fa.gz"
                        handle.write(
                            "mv "                                        + \
                            outdir+"/query.edit.gz."+str(i)+".gz.fa.gz " + \
                            moved + "." + str(i) + ".gz\n"
                            )


                    if(ML_or_MP=="ML"): 
                        if(ALIGNED=="unaligned"): # for unaligned sequences
                            # Conduct HMM alignment
                            handle.write(
                                hmm_aligner                                +
                                " --outformat afa"                         +
                                #" -q "                                     +
                                #" -m "                                     +
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
                                " -in " + outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa "          +
                                "-selectcols { "                                           +
                                "   `trimal -sgc "                                         +
                                "    -in " +outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa.ref "    +
                                "    |awk ' { if( $2==100 ){ print $1 }}'"                 +
                                "    |tr \"\\n\" \",\" ` "                                 +
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
                            seed                                         +
                            "\n"
                            )
                    elif(ML_or_MP=="MP"):
                        handle.write(
                            "cd "+outdir+"/EPANG"+str(i)+"\n"
                        )
                        if(ALIGNED=="unaligned"): # for unaligned sequences
                            # Conduct HMM alignment
                            handle.write(
                                hmm_aligner                                 +
                                " --outformat afa"                          +
                                " --mapali "+refseq+" "                     +
                                refseq+".hmm "+moved+"."+str(i)+".gz"       +
                                " | sed 's/\./N/g'"                         +
                                " > "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa\n"
                            )   
                            handle.write(
                                "rm " + queryfile + "\n"
                            )   
                        elif(ALIGNED=="aligned"): # for aligned sequences
                            handle.write(
                                "cat "+refseq+" "+
                                moved+"."+str(i)+".gz"+
                                " | gunzip > "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa\n"
                            )
                            handle.write(
                                "rm " + moved+"."+str(i)+".gz" + "\n"
                            )  
                        handle.write(
                            RAXMLSEQ                                      +
                            " -n epa_result -f y -m GTRCAT"               +
                            " -s "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa" +
                            " -t "+reftree+"\n"
                            ) 
                        handle.write(
                            "rm "+outdir+"/EPANG"+str(i)+"/"+filename+"/ref_query.fa*\n"
                            ) 
                        handle.write(
                            "python3 "                                                      +
                            codedir+"/python/jplace_parse.py "                              +
                            outdir+"/EPANG"+str(i)+"/"+filename+"/RAxML_portableTree.epa_result.jplace " +
                            "epa_MP "                                                       +
                            seed + "\n"
                            )
                handle.write(
                    "echo \"finished\" > "      +
                    outdir+"/epang"+str(i)+".o\n"
                    )
                #handle.write(
                #    "rm "+outdir+"/*."+str(i)+".gz\n"
                #    ) 
                # end of a distributed task
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
        
        if (file_format == "fasta"):
            shutil.move(moved,query)
        elif (file_format == "edit"):
            shutil.move(outdir+"/query.edit.gz",query)

        # (If HMM alignments were conducted) concat aligned sequences
        if(ALIGNED=="unaligned"):
            subprocess.call(
                "cat "+
                outdir+"/EPANG*/*/ref_query.fa.selectcols.query "+
                " | gzip > "+
                WD+"/INPUT.fa.aligned.gz",
                shell=True
            )
            subprocess.call(
                "gzip -c "+
                "$(ls "+outdir+"/EPANG0/*/ref_query.fa.selectcols.ref | head -n1)"+
                "> "+ 
                WD+"/SUBSAMPLE.fa.aligned.gz",
                shell=True
                )

        # merge results
        #shutil.move(outdir+"/EPANG0/placement_tree.out",outdir+"/placement_tree.out")
        #my_paste(outdir,nodenum, outdir+"/edge_to_seqname_all.out")
        subprocess.call(
            "mv "  + "$(ls "+outdir+"/EPANG0/*/placement_tree.out | head -n1) "+outdir+"/placement_tree.out;"+
            "cat " + outdir + "/EPANG*/*/edge_to_seqname.out > " + outdir+"/edge_to_seqname_all.out",
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