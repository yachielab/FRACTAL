echo -n "Lineage reconstruction started:  "
date

if [ $# -ne 18 ]; then
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
SEED=${17}
JOB_NAME=${18}
ROOT_DIR=${DATA_DIR}/${exp_num}
QSUB_OPTION="${QSUB_OPTION} -N ${JOB_NAME}"

# QSUB OPTION
if [ -z "$INIT_QSUB_OPTION" ]; then
  INIT_QSUB_OPTION=$QSUB_OPTION
fi

mkdir ${ROOT_DIR}
mkdir ${ROOT_DIR}/nodes
mkdir ${ROOT_DIR}/qsub_dir
mkdir ${ROOT_DIR}/final_tree
mkdir ${ROOT_DIR}/out
mkdir ${ROOT_DIR}/err
mkdir ${ROOT_DIR}/executed

# requirement check
EPANG=$(which epa-ng)
RAXMLSEQ=$(which raxmlHPC-SSE3)
RAXMLPAR=$(which raxmlHPC-PTHREADS-SSE3)

#setting for the 1st qsub
mkdir ${ROOT_DIR}/nodes/d0
cp ${input_faname} ${ROOT_DIR}/nodes/d0/INPUT.fa
while [ ! -e ${ROOT_DIR}/nodes/d0/INPUT.fa ]; do
  echo "${ROOT_DIR}/nodes/d0/INPUT.fa was not found" 1>&2
  cp ${input_faname} ${ROOT_DIR}/nodes/d0/INPUT.fa
done

echo "1" >${ROOT_DIR}/NUMFILE
echo "#!/bin/bash" >${ROOT_DIR}/qsub_dir/qsub_d0.sh
echo "#$ -S /bin/bash" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
echo "export PATH=${PATH}" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh
echo "python3 ${CODE_DIR}/python/FRACluster.py ${ROOT_DIR}/nodes/d0 ${num_of_subsample} ${subsample_size} ${ROOT_DIR}/nodes $threshold ${THREADNUM} ${ROOT_DIR}/NUMFILE ${ROOT_DIR}/qsub_dir ${CODE_DIR} $ROOTING $MODEL \"${OPTION}\" ${TREE} aligned $EPANG $RAXMLSEQ $RAXMLPAR $SOFTWARE $max_num_of_jobs 0 \"$SEED\"" >>${ROOT_DIR}/qsub_dir/qsub_d0.sh

# first qsub
if [ $max_num_of_jobs -gt 1 ]; then
  qsub ${INIT_QSUB_OPTION} -o ${ROOT_DIR}/out -e ${ROOT_DIR}/err ${ROOT_DIR}/qsub_dir/qsub_d0.sh # parallel mode
  wait
else
  bash ${ROOT_DIR}/qsub_dir/qsub_d0.sh >${ROOT_DIR}/out/qsub_d0.sh.out 2>${ROOT_DIR}/err/qsub_d0.sh.err # sequential mode
  wait
fi
mv ${ROOT_DIR}/qsub_dir/qsub_d0.sh ${ROOT_DIR}/executed/qsub_d0.sh
sleep 5

# keep looking over qsub directory & executing qsub until no job is being processed in calculation node

if [ $max_num_of_jobs -gt 1 ]; then # parallel mode
  a=$(qstat | grep ${JOB_NAME} | wc -l)
  b=$(ls ${ROOT_DIR}/qsub_dir | wc -l)
  while [ $(expr $a + $b) -gt 0 ]; do
    for file in $(ls ${ROOT_DIR}/qsub_dir); do
      NUMBER_OF_JOBS=$(qstat | grep ${JOB_NAME} | wc -l)
      wait
      if [ -z $NUMBER_OF_JOBS ]; then NUMBER_OF_JOBS=0; fi
      if [ $NUMBER_OF_JOBS -lt ${max_num_of_jobs} ]; then
        qsub ${QSUB_OPTION} -o ${ROOT_DIR}/out -e ${ROOT_DIR}/err ${ROOT_DIR}/qsub_dir/${file}
        wait
        mv ${ROOT_DIR}/qsub_dir/${file} ${ROOT_DIR}/executed/${file}
      fi
    done
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
  qsub ${QSUB_OPTION} -o ${ROOT_DIR}/out -e ${ROOT_DIR}/err ${ROOT_DIR}/qsub_dir/qsub_assembly.sh
  wait
else
  bash ${ROOT_DIR}/qsub_dir/qsub_assembly.sh
  wait
fi
mv ${ROOT_DIR}/qsub_dir/qsub_assembly.sh ${ROOT_DIR}/executed/qsub_assembly.sh

while [ ! -e ${ROOT_DIR}/final_tree/assembly_flag.txt ]; do
  sleep 5
done

echo -n "Lineage reconstruction finished: "
date
