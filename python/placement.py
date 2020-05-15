from Bio import SeqIO
import jplace_parse
import divide_ref_and_query
import subprocess
import os
import shutil

# decompose fasta file (in_file) into x fasta files.
def decompose_fasta(in_file, x,seq_count):
    n=seq_count
    k=n//x + 1 # each decomposed file has k sequences
    ohandle=[]
    for i in range(x):
        ohandle.append(
            open(in_file+"."+str(i), 'w')
        ) 
    with open(in_file,'r') as ihandle:
        allseq_itr = SeqIO.parse(ihandle, "fasta")
        l=0 # inclement constantly
        m=0 # inclement for each sequence, but become 0 after m reaches k
        i=0 # 
        while(l<n):
            for record in allseq_itr:
                if(m<k):
                    m+=1
                else:
                    i+=1
                    m=0
                SeqIO.write(record,ohandle[i],'fasta')
                l+=1
    for i in range(x):
        ohandle[i].close()

def distributed_placement(  WD, EPANG, refseq, reftree, model, 
                            query, outdir, threadnum, nodenum, 
                            codedir, seq_count, ML_or_MP, RAXMLSEQ, 
                            ALIGNED, seed, hmm_aligner="", hmm_profiler=""  ):
    if(nodenum<=1):
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
                    "| sed 's/\./N/g'> "   +
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
                shutil.move(outdir+"/ref_query.fa.query", WD+"/INPUT.fa.aligned")
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
                seed
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
                    outdir+"/ref_query.fa"     ,
                    shell=True
                )  
            elif(ALIGNED=="aligned"): # for aligned sequences
                subprocess.call(
                    "cat "+refseq+" "            +
                    query                        +
                    " > "+outdir+"/ref_query.fa" ,
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
                seed
            )
        os.rename(
            outdir+"/edge_to_seqname.out",
            outdir+"/edge_to_seqname_all.out"
            )
    else:
        dname=WD.split("/").pop()
        moved=outdir+"/query.fa"
        shutil.move(query, moved)
        # Build HMM profile
        subprocess.call(
            hmm_profiler+" "+
            refseq+".hmm "+
            refseq,
            shell=True
            )
        decompose_fasta(
            moved,
            nodenum, 
            seq_count
            )

        #distribution start
        for i in range(nodenum):
            os.mkdir(outdir+"/EPANG"+str(i))
            with open(WD+"/../../qsub_dir/qsub_"+dname+"."+str(i)+".sh", 'w') as handle:
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
                if(ML_or_MP=="ML"): 
                    if(ALIGNED=="unaligned"): # for unaligned sequences
                        # Conduct HMM alignment
                        handle.write(hmm_aligner+
                            " --outformat Phylip"                      +
                            " -m "                                     +
                            " --trim "                                 +
                            " --mapali " + refseq + " "                +
                            refseq+".hmm "                             +
                            moved +"."+str(i)                          +
                            " | sed 's/\./N/g'> "                      +
                            outdir+"/EPANG"+str(i)+"/ref_query.fa\n"
                            )   
                        handle.write(
                            "python3 "                                 +
                            codedir+"/python/divide_ref_and_query.py " +
                            outdir+"/EPANG"+str(i)+"/ref_query.fa "    + 
                            refseq + "\n"
                            )
                        handle.write(
                            EPANG                                               + 
                            " --redo"                                           +
                            " -s "+outdir+"/EPANG"+str(i)+"/ref_query.fa.ref"   +
                            " -t "+reftree                                      +
                            " --model "+model                                   +
                            " -q "+outdir+"/EPANG"+str(i)+"/ref_query.fa.query" +
                            " -w "+outdir+"/EPANG"+str(i)                       +
                            " -T "+str(threadnum)+"\n"
                            )
                    elif(ALIGNED=="aligned"): # for aligned sequences
                        handle.write(
                            EPANG                         +
                            " --redo"                     +
                            " -s "+refseq                 +
                            " -t "+reftree                +
                            " --model "+model             +
                            " -q "+moved+"."+str(i)       +
                            " -w "+outdir+"/EPANG"+str(i) +
                            " -T "+str(threadnum)+"\n"
                            )
                    handle.write(
                        "cd "+outdir+"/EPANG"+str(i)+"\n"
                        )
                    handle.write(
                        "python3 "                                   +
                        codedir+"/python/jplace_parse.py "           +
                        outdir+"/EPANG"+str(i)+"/epa_result.jplace " +
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
                            hmm_aligner+
                            " --outformat afa"+
                            " --mapali "+refseq+" "+
                            refseq+".hmm "+moved+"."+str(i)+
                            " | sed 's/\./N/g'> "+
                            outdir+"/EPANG"+str(i)+"/ref_query.fa\n"
                        )   
                    elif(ALIGNED=="aligned"): # for aligned sequences
                        handle.write(
                            "cat "+refseq+" "+
                            moved+"."+str(i)+
                            " > "+outdir+"/EPANG"+str(i)+"/ref_query.fa\n"
                        )
                    handle.write(
                        RAXMLSEQ                                      +
                        " -n epa_result -f y -m GTRCAT"               +
                        " -s "+outdir+"/EPANG"+str(i)+"/ref_query.fa" +
                        " -t "+reftree+"\n"
                        ) 
                    handle.write(
                        "python3 "                                                      +
                        codedir+"/python/jplace_parse.py "                              +
                        outdir+"/EPANG"+str(i)+"/RAxML_portableTree.epa_result.jplace " +
                        "epa_MP "                                                       +
                        seed + "\n"
                        )
                handle.write(
                    "echo \"finished\" > "      +
                    outdir+"/epang"+str(i)+".o"
                    )
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
            os.remove(outdir+"/query.fa."+str(i))
        
        shutil.move(moved,query)

        # (If HMM alignments were conducted) concat aligned sequences
        if(ALIGNED=="unaligned"):
            subprocess.call(
                "cat "+
                outdir+"/EPANG*/ref_query.fa.query "+
                "> "+
                WD+"/INPUT.fa.aligned",
                shell=True
            )

        # merge results
        shutil.move(outdir+"/EPANG0/placement_tree.out",outdir+"/placement_tree.out")
        my_paste(outdir,nodenum, outdir+"/edge_to_seqname_all.out")

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