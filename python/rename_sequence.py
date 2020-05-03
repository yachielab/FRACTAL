from Bio import SeqIO
import random
import sys
import shutil
import os
import subprocess
import jplace_parse

def rename_sequence(in_fname,out_fname):
    with open(in_fname) as origin, open(out_fname, 'w') as renamed:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k = 0
        for s in input_itr:
            s.id = "s"+str(k)
            s.description = "s"+str(k)
            k = k +1
            SeqIO.write(s, renamed, "fasta")


def outgroup_check_fast(in_fname):
    exist_root = False
    with open(in_fname) as origin:
        for line in origin:
            if (line==">root\n"): exist_root = True; break
    return exist_root

def count_sequence(in_fname):
    with open(in_fname) as origin:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        k=0
        for s in input_itr:
            k+=1
        input_itr = SeqIO.parse(origin, "fasta")
        n=0
        for s in input_itr:
            n+=len(str(s.seq))//k
    return [k,n] # k: number of sequence, n: average sequence length

def count_sequence_fast(in_fname):
    with open(in_fname) as handle:
        k, l=0, 0
        for line in handle:
            if(line[0]==">"): k+=1
            elif(k==1): l+=len(line)-1
    return [k,l] # k: number of sequence, n: sequence length of first sequence (outgroup)

def random_sampling(in_fname,out_fname,subsample_size,seed,n=None):
    if(seed=="random"): random.seed(int(random.randint(0,99999)))
    elif(len(seed)!=0):random.seed(int(seed))
    else:print("-r Error: invalid random seed!")
    if (n==None):
        n=count_sequence_fast(in_fname)[0]
    rand_idx=random.sample(range(n-1),subsample_size)
    rand_idx.sort()
    with open(out_fname, 'w') as subs:
        with open(in_fname) as allseq:
            allseq_itr = SeqIO.parse(allseq, "fasta")
            for s in allseq_itr:
                if(s.id=="root"):
                    SeqIO.write(s, subs, "fasta")
        with open(in_fname) as allseq:
            allseq_itr = SeqIO.parse(allseq, "fasta")
            i=0 # index on rand_idx
            k=0 # index on record
            for s in allseq_itr:
                if(i>=len(rand_idx)):
                    break
                if(s.id!="root"):
                    if(k==rand_idx[i]):
                        SeqIO.write(s, subs, "fasta")
                        i+=1
                    k+=1

# decompose fasta file (in_file) into x fasta files.
def decompose_fasta(in_file, x,seq_count):
    n=seq_count
    k=n//x + 1 # each decomposed file has k sequences
    ohandle=[]
    for i in range(x):
        ohandle.append(open(in_file+"."+str(i), 'w')) 
    with open(in_file,'r') as ihandle:
        allseq_itr = SeqIO.parse(ihandle, "fasta")
        l=0 # inclement constantly
        m=0 # inclement for each sequence, but become 0 after m reaches k
        i=0 # 
        while(l<n):
            sys.stdout.write("\rsequence: %d" % l)
            sys.stdout.flush()
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

def distributed_placement(WD, EPANG, refseq, reftree, model, query, outdir, threadnum, nodenum, codedir, seq_count, ML_or_MP, RAXMLSEQ):
    if(nodenum<=1):
        if(ML_or_MP=="ML"): 
            subprocess.call(EPANG+" --redo -s "+refseq+" -t "+reftree+" --model "+model+" -q "+query+" -w "+outdir+" -T "+str(threadnum),shell=True)
            os.chdir(outdir)
            jplace_parse.parse_jplace(outdir+"/epa_result.jplace",placement_method="epa-ng")
        if(ML_or_MP=="MP"): 
            subprocess.call("cat "+refseq+" "+query+" > "+outdir+"/ref_query.fa",shell=True)
            os.chdir(outdir)
            subprocess.call(RAXMLSEQ+" -n epa_result -f y -m GTRCAT -s "+outdir+"/ref_query.fa"+" -t "+reftree,shell=True)
            jplace_parse.parse_jplace(outdir+"/RAxML_portableTree.epa_result.jplace",placement_method="epa_MP")
        os.rename(outdir+"/edge_to_seqname.out", outdir+"/edge_to_seqname_all.out")
    else:
        dname=WD.split("/").pop()
        moved=outdir+"/query.fa"
        shutil.move(query, moved)
        decompose_fasta(moved, nodenum, seq_count)

        #distribution start
        for i in range(nodenum):
            os.mkdir(outdir+"/EPANG"+str(i))
            with open(WD+"/../../qsub_dir/qsub_"+dname+"."+str(i)+".sh", 'w') as handle:
                PATH = (subprocess.Popen('echo $PATH', stdout=subprocess.PIPE,shell=True).communicate()[0]).decode('utf-8')
                PATH = (PATH.split('\n'))[0]
                LD_LIBRARY_PATH = (subprocess.Popen('echo $LD_LIBRARY_PATH', stdout=subprocess.PIPE,shell=True).communicate()[0]).decode('utf-8')
                LD_LIBRARY_PATH = (LD_LIBRARY_PATH.split('\n'))[0]
                handle.write("#!/bin/bash\n")
                handle.write("#$ -S /bin/bash\n")
                handle.write("PATH={}\n".format(PATH))
                handle.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
                if(ML_or_MP=="ML"): 
                    handle.write(EPANG+" --redo -s "+refseq +" -t "+reftree+" --model "+model+" -q "+moved+"."+str(i)+" -w "+outdir+"/EPANG"+str(i)+" -T "+str(threadnum)+"\n")
                    handle.write("cd "+outdir+"/EPANG"+str(i)+"\n")
                    handle.write("python3 "+codedir+"/python/jplace_parse.py "+outdir+"/EPANG"+str(i)+"/epa_result.jplace epa-ng\n")
                elif(ML_or_MP=="MP"):
                    subprocess.call("cat "+refseq+" "+query+" > "+outdir+"/ref_query.fa",shell=True)
                    handle.write("cd "+outdir+"/EPANG"+str(i)+"\n")
                    subprocess.call(RAXMLSEQ+" -n epa_result -f y -m GTRCAT -s "+outdir+"/ref_query.fa"+" -t "+reftree,shell=True) 
                    handle.write("python3 "+codedir+"/python/jplace_parse.py "+outdir+"/EPANG"+str(i)+"/epa_result.jplace epa_MP\n")
                handle.write("echo \"finished\" > "+outdir+"/epang"+str(i)+".o")
        #distribution end
        flag=0
        while(flag==0):
            i=0
            while(i<nodenum):
                if(not os.path.exists(outdir+"/epang"+str(i)+".o")):
                    break
                i+=1
            if i == nodenum:
                flag=1
        for i in range(nodenum):
            os.remove(outdir+"/epang"+str(i)+".o")
            os.remove(outdir+"/query.fa."+str(i))
        shutil.move(moved,query)
        shutil.move(outdir+"/EPANG0/placement_tree.out",outdir+"/placement_tree.out")
        #subprocess.call("paste "+outdir+"/EPANG*/edge_to_seqname.out | tr -d '\t'> "+outdir+"/edge_to_seqname_all.out",shell=True)
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

'''
command line argument: "<input .fa file path> <output .fa file path>"
'''
if __name__ == "__main__":
    argvs = sys.argv
    rename_sequence(argvs[1],argvs[2])