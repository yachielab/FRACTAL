CMDNAME="TREE.sh"

THREAD_NUM=4
TREEMETHOD="unspecified"
ALIGNED='aligned'
FILE='unspecified'
CODE_DIR='unspecified'
WD='unspecified' # working directory
OPTION=""
MODEL="GTRCAT" # model of site heterogeneity in raxml
SOFTWARE="unspecified"
ALIGNER="unspecified"

while getopts  n:m:a:f:c:w:p:d:q:s: OPT
do
  case $OPT in
    "n" ) FLG_TN="TRUE" ; THREADNUM="$OPTARG";;
    "m" ) FLG_TM="TRUE" ; TREEMETHOD="$OPTARG" ;;
    "a" ) FLG_AM="TRUE" ; ALIGNED="$OPTARG" ;;
    "f" ) FLG_F="TRUE" ; FILE="$OPTARG" ;;
    "c" ) FLG_C="TRUE" ; CODE_DIR="$OPTARG" ;;
    "w" ) FLG_WD="TRUE" ; WD="$OPTARG" ;;
    "p" ) FLG_OPT="TRUE" ; OPTION="$OPTARG" ;;
    "d" ) FLG_MOD="TRUE" ; MODEL="$OPTARG" ;;
    "q" ) FLG_NJ="TRUE" ; SOFTWARE="$OPTARG" ;;
    "s" ) FLG_S="TRUE" ; ALIGNER="$OPTARG" ;;
      * ) echo "Usage:" 1>&2
          echo $CMDNAME 1>&2
          echo "[-n number of threads]" 1>&2
          echo "[-m tree estimation method]" 1>&2
          echo "[-a aligned or not]" 1>&2
          echo "[-f input .fa file path]" 1>&2
          echo "[-c code directory]" 1>&2
          echo "[-w working directory]" 1>&2
          echo "[-p option]" 1>&2
          echo "[-d evolution model of ML]" 1>&2
          echo "[-q tree construction software executable path]" 1>&2
          exit 1 ;;
  esac
done

if [ "$OPTION" = "nothing" ]; then
    OPTION=""
fi

cd $WD

# Alignment
if [ $ALIGNED = 'aligned' ]; then
    echo "Already aligned!"
    gunzip -c ${FILE} > ${FILE}.aligned
elif [ $ALIGNED = "unaligned" ]; then
    ${ALIGNER} --quiet <(gunzip -c ${FILE}) > ${FILE}.aligned # set $PATH!!!!!!!!!!
else
    echo "exception: Alignment method name seems wrong..."
fi

# Tree construction
# ML (RAxML)
if [ $TREEMETHOD = "raxmlML" ]; then
    ${SOFTWARE} -s ${FILE}.aligned -n raxml -T $THREAD_NUM ${OPTION} -m ${MODEL} -p 12345
    wait
    mv RAxML_bestTree.raxml ${FILE}.aligned.tree
# NJ
elif [ $TREEMETHOD = "rapidnjNJ" ]; then
    $SOFTWARE -c ${THREADNUM} ${FILE}.aligned -i fa ${OPTION} > ${FILE}.aligned.tree
# MP
elif [ $TREEMETHOD = "raxmlMP" ]; then
    ${SOFTWARE} -s ${FILE}.aligned -y -n raxml -T $THREAD_NUM -m ${MODEL} -p 12345 # -y : only compute maximum parsimony!
    mv RAxML_parsimonyTree.raxml ${FILE}.aligned.tree
# ML (fasttree)
elif [ $TREEMETHOD = "fasttreeML" ]; then
    export OMP_NUM_THREADS=$THREADNUM
    $SOFTWARE -gtr -nt -seed 111 -quiet ${OPTION}  < ${FILE}.aligned > ${FILE}.aligned.tree
# Others
elif [ $TREEMETHOD = "unspecified" ]; then
    $SOFTWARE ${FILE}.aligned
else
    echo "exception: Tree construction method name seems wrong..."
fi

#gzip ${FILE}.aligned