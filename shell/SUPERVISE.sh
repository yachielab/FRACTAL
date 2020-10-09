echo -n "Lineage reconstruction started:  "
date

if [ $# -ne 27 ]; then
  echo "args:$#" 1>&2
  echo "SUPERVISE.sh: wrong number of arguments!" 1>&2
  exit 1
fi

num_of_subsample=$1 #Integer
subsample_size=$2   #Integer
input_faname=$3     #Absolute Path of .fa file
exp_num=$4          #example: "NK_0001"
CODE_DIR=$5         #"......./code"
DATA_DIR=$6         #"......./data"
threshold=$7
ROOTING=$8
max_num_of_jobs=$9
TREE=${10}
THREADNUM=${11}
SOFTWARE=${12}
OPTION=${13}
MODEL=${14}
QSUB_OPTION=${15}
INIT_QSUB_OPTION=${16}
ASSEMBLY_QSUB_OPTION=${17}
SEED=${18}
JOB_NAME=${19}
PLACEMENT_METHOD=${20}
ALIGNED=${21}
max_num_of_iterations=${22}
extraction_size=${23}
careful=${24}
FASTA_or_EDIT=${25}
GZIP_INTERMEDIATE=${26}
SEQ_NUM_FILE=${27}
ROOT_DIR=${DATA_DIR}/${exp_num}


# QSUB OPTION
if [ -z "$INIT_QSUB_OPTION" ]; then
  INIT_QSUB_OPTION=$QSUB_OPTION
fi
if [ -z "$ASSEMBLY_QSUB_OPTION" ]; then
  ASSEMBLY_QSUB_OPTION=$QSUB_OPTION
fi

mkdir ${ROOT_DIR}
mkdir ${ROOT_DIR}/nodes
mkdir ${ROOT_DIR}/qsub_dir
mkdir ${ROOT_DIR}/final_tree
mkdir ${ROOT_DIR}/out
mkdir ${ROOT_DIR}/err
mkdir ${ROOT_DIR}/executing
mkdir ${ROOT_DIR}/executed

# requirement check
EPANG=$(which epa-ng)
RAXMLSEQ=$(which raxmlHPC-SSE3)
RAXMLPAR=$(which raxmlHPC-PTHREADS-SSE3)

if [ "$ALIGNED" = "unaligned" ]; then
    MAFFT=$(which mafft)
    HMM_BUILD=$(which hmmbuild)
    HMM_ALIGN=$(which hmmalign)
fi

#<<COMMENTOUT

# setting for the 1st qsub
mkdir ${ROOT_DIR}/nodes/d0

if   [ -d ${input_faname} ]; then
    filepath_list=$(ls ${input_faname}/*)
elif [ -f ${input_faname} ]; then
    filepath_list=${input_faname}
fi

for input_faname in ${filepath_list}; do 
    # input gzipped or not
    if [ $(echo ${input_faname} | sed 's/^.*\.\([^\.]*\)$/\1/') = "gz" ]; then
        gzip_input="gunzip"
    else
        gzip_input="cat"
    fi

    # output gzipped or not
    if [ $GZIP_INTERMEDIATE = "TRUE" ]; then
        gzip_output="gzip"
        out_extention=".gz"
    else
        gzip_output="cat"
        out_extention=""
    fi

    # avoid gunzip & gzip
    if [ $gzip_input = "gunzip" -a $gzip_output = "gzip" ]; then
        gzip_input="cat"
        gzip_output="cat"
    fi

    copied_fpath=${ROOT_DIR}/nodes/d0/$(basename ${input_faname}).${FASTA_or_EDIT}${out_extention}

    cat ${input_faname} | ${gzip_input} | ${gzip_output} > $copied_fpath

    if [ -e $SEQ_NUM_FILE ]; then
        (echo -ne "${copied_fpath}\t"; cat $SEQ_NUM_FILE | grep $(basename ${input_faname}) | cut -f2) >> ${ROOT_DIR}/nodes/d0/file2Nseq.txt
    fi

done

wait
echo "1" >${ROOT_DIR}/NUMFILE
echo "#!/bin/bash" >${ROOT_DIR}/qsub_dir/qsub_d0.sh
echo "#$ -S /bin/bash" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
echo "export PATH=${PATH}" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
if [ "$FASTA_or_EDIT" = "fa" ]; then
    echo "python3 ${CODE_DIR}/python/FRACluster.py   ${ROOT_DIR}/nodes/d0 ${num_of_subsample} ${subsample_size} ${ROOT_DIR}/nodes $threshold ${THREADNUM} ${ROOT_DIR}/NUMFILE ${ROOT_DIR}/qsub_dir ${CODE_DIR} $ROOTING $MODEL \"${OPTION}\" ${TREE} ${ALIGNED} $EPANG $RAXMLSEQ $RAXMLPAR $SOFTWARE $max_num_of_jobs 0 \"$SEED\" ${PLACEMENT_METHOD} ${extraction_size} ${careful} ${MAFFT} ${HMM_BUILD} ${HMM_ALIGN} 0" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
elif [ "$FASTA_or_EDIT" = "edit" ]; then
    echo "python3 ${CODE_DIR}/python/FRACluster.2.py ${ROOT_DIR}/nodes/d0 ${num_of_subsample} ${subsample_size} ${ROOT_DIR}/nodes $threshold ${THREADNUM} ${ROOT_DIR}/NUMFILE ${ROOT_DIR}/qsub_dir ${CODE_DIR} $ROOTING $MODEL \"${OPTION}\" ${TREE} ${ALIGNED} $EPANG $RAXMLSEQ $RAXMLPAR $SOFTWARE $max_num_of_jobs 0 \"$SEED\" ${PLACEMENT_METHOD} ${extraction_size} ${careful} ${MAFFT} ${HMM_BUILD} ${HMM_ALIGN} 0" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
fi

# first qsub
if [ $max_num_of_jobs -gt 1 ]; then
  qsub_err="yet"
  while [ -n "${qsub_err}" ]; do
    if [ $qsub_err != "yet" ]; then
        echo ${qsub_err}
    fi
    qsub_err=$((qsub -N ${JOB_NAME} ${INIT_QSUB_OPTION} -o ${ROOT_DIR}/out/qsub_d0.sh.out -e ${ROOT_DIR}/err/qsub_d0.sh.err ${ROOT_DIR}/qsub_dir/qsub_d0.sh 1> /dev/null) 2>&1) # parallel mode
    wait
  done
else
  bash ${ROOT_DIR}/qsub_dir/qsub_d0.sh >${ROOT_DIR}/out/qsub_d0.sh.out 2>${ROOT_DIR}/err/qsub_d0.sh.err # sequential mode
  wait
fi
echo "qsub... ${ROOT_DIR}/qsub_dir/qsub_d0.sh submitted!" 
mv ${ROOT_DIR}/qsub_dir/qsub_d0.sh ${ROOT_DIR}/executed/qsub_d0.sh
sleep 5

#COMMENTOUT

# keep looking over qsub directory & executing qsub until no job is being processed in calculation node

if [ $max_num_of_jobs -gt 1 ]; then # parallel mode
  a=$(qstat | grep ${JOB_NAME} | wc -l)
  b=$(ls ${ROOT_DIR}/qsub_dir | wc -l)
  while [ $(expr $a + $b) -gt 0 ]; do
    for file in $(ls ${ROOT_DIR}/qsub_dir/*placement* 2>/dev/null) $(ls ${ROOT_DIR}/qsub_dir/*partition* 2>/dev/null) $(ls ${ROOT_DIR}/qsub_dir/*cycle* 2>/dev/null); do
      NUMBER_OF_JOBS=$(qstat | grep ${JOB_NAME} | wc -l)
      wait
      if [ -z $NUMBER_OF_JOBS ]; then NUMBER_OF_JOBS=0; fi
      if [ $NUMBER_OF_JOBS -lt ${max_num_of_jobs} ]; then
        qsub_err="yet"
        while [ -n "${qsub_err}" ]; do
            if [ "$qsub_err" != "yet" ]; then
                echo ${qsub_err}
            fi
            if [ `echo ${file} | grep 'largemem'` ] ; then
                qsub_err=$((qsub ${INIT_QSUB_OPTION} -N ${JOB_NAME} -o ${ROOT_DIR}/out/${file}.out -e ${ROOT_DIR}/err/${file}.err ${ROOT_DIR}/qsub_dir/${file} 1> /dev/null) 2>&1)
            else
                qsub_err=$((qsub ${QSUB_OPTION}      -N ${JOB_NAME} -o ${ROOT_DIR}/out/${file}.out -e ${ROOT_DIR}/err/${file}.err ${ROOT_DIR}/qsub_dir/${file} 1> /dev/null) 2>&1)
            fi
        done
        wait
        echo "qsub... ${ROOT_DIR}/qsub_dir/${file} submitted!" 
        mv ${ROOT_DIR}/qsub_dir/${file} ${ROOT_DIR}/executed/${file}
      fi
    done
    NUMBER_OF_ITERATIONS=$(ls ${ROOT_DIR}/nodes | wc -l)

    # for safety
    if [ $NUMBER_OF_ITERATIONS -gt ${max_num_of_iterations} ]; then
        echo "Number of FRACTAL iterations > "${max_num_of_iterations}"!"
        echo "FRACTAL gives up the lineage reconstruction to avoid generating too much intermediate files."
        echo "If you would like to raise the upper limitation of FRACTAL iterations,"
        echo "  please set the limit by -l option, then please carefully supervise"
        echo "  the number of intermediate files generated by FRACTAL during calculation."
        exit 1
    fi
    
    sleep 5
    a=$(qstat | grep ${JOB_NAME} | wc -l)
    b=$(ls ${ROOT_DIR}/qsub_dir | wc -l)
  done
else # sequential mode
  while [ $(ls ${ROOT_DIR}/qsub_dir | wc -l) -gt 0 ]; do
    for file in $(ls ${ROOT_DIR}/qsub_dir); do
      bash ${ROOT_DIR}/qsub_dir/${file} >${ROOT_DIR}/out/${file}.out 2>${ROOT_DIR}/err/${file}.err
      wait
      mv ${ROOT_DIR}/qsub_dir/${file} ${ROOT_DIR}/executed/${file}
    done
    NUMBER_OF_ITERATIONS=$(ls ${ROOT_DIR}/nodes | wc -l)

    # for safety
    if [ $NUMBER_OF_ITERATIONS -gt ${max_num_of_iterations} ]; then
        echo "Number of FRACTAL iterations > "${max_num_of_iterations}"!"
        echo "FRACTAL gives up the lineage reconstruction to avoid generating too much intermediate files."
        echo "If you would like to raise the upper limitation of FRACTAL iterations,"
        echo "  please set the limit by -l option, then please carefully supervise"
        echo "  the number of intermediate files generated by FRACTAL during calculation."
        exit 1
    fi

    sleep 5
  done
fi

echo -n "All FRACTAL iterations finished: "
date

# tree infetence ended here

# assemble the subtree .nwk files
echo "#!/bin/bash" >${ROOT_DIR}/qsub_dir/qsub_assembly.sh
echo "#$ -S /bin/bash" >>${ROOT_DIR}/qsub_dir/qsub_assembly.sh
echo "export PATH=${PATH}" >>${ROOT_DIR}/qsub_dir/qsub_assembly.sh
echo "python3 ${CODE_DIR}/python/TreeAssembly.py ${ROOT_DIR}/nodes/d0 ${ROOT_DIR}/final_tree/HUGE_Result.nwk TRUE" >>${ROOT_DIR}/qsub_dir/qsub_assembly.sh
echo "echo \"finished\" > ${ROOT_DIR}/final_tree/assembly_flag.txt" >>${ROOT_DIR}/qsub_dir/qsub_assembly.sh
if [ $max_num_of_jobs -gt 1 ]; then
    qsub_err="yet"
    while [ -n "${qsub_err}" ]; do
        qsub_err=$((qsub ${ASSEMBLY_QSUB_OPTION} -N ${JOB_NAME} -o ${ROOT_DIR}/out/qsub_assembly.sh.out -e ${ROOT_DIR}/err/qsub_assembly.sh.err ${ROOT_DIR}/qsub_dir/qsub_assembly.sh 1> /dev/null) 2>&1)
        wait
    done
else
  bash ${ROOT_DIR}/qsub_dir/qsub_assembly.sh
  wait
fi
echo "qsub... ${ROOT_DIR}/qsub_dir/qsub_assembly.sh submitted!" 
mv ${ROOT_DIR}/qsub_dir/qsub_assembly.sh ${ROOT_DIR}/executed/qsub_assembly.sh

while [ ! -e ${ROOT_DIR}/final_tree/assembly_flag.txt ]; do
  sleep 5
done

echo -n "Lineage reconstruction finished: "
date