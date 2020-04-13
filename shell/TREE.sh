CMDNAME="TREE.sh"

THREAD_NUM=4
TREEMETHOD="unspecified"
ALIGNMETHOD='aligned'
FILE='unspecified'
CODE_DIR='unspecified'
WD='unspecified' # working directory
OPTION=""
MODEL="GTRCAT" # model of site heterogeneity in raxml
SOFTWARE="unspecified"

while getopts  n:m:a:f:c:w:p:d:q:r: OPT
do
  case $OPT in
    "n" ) FLG_TN="TRUE" ; THREADNUM="$OPTARG";;
    "m" ) FLG_TM="TRUE" ; TREEMETHOD="$OPTARG" ;;
    "a" ) FLG_AM="TRUE" ; ALIGNMETHOD="$OPTARG" ;;
    "f" ) FLG_F="TRUE" ; FILE="$OPTARG" ;;
    "c" ) FLG_C="TRUE" ; CODE_DIR="$OPTARG" ;;
    "w" ) FLG_WD="TRUE" ; WD="$OPTARG" ;;
    "p" ) FLG_OPT="TRUE" ; OPTION="$OPTARG" ;;
    "d" ) FLG_MOD="TRUE" ; MODEL="$OPTARG" ;;
    "q" ) FLG_NJ="TRUE" ; SOFTWARE="$OPTARG" ;;
      * ) echo "Usage:" 1>&2
          echo $CMDNAME 1>&2
          echo "[-n number of threads]" 1>&2
          echo "[-m tree estimation method]" 1>&2
          echo "[-a alignment method]" 1>&2
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
echo $ALIGNMETHOD
if [ $ALIGNMETHOD = 'aligned' ]; then
    echo "Already aligned!"
    mv ${FILE} ${FILE}.aligned
elif [ $ALIGNMETHOD = "mafft" ]; then
    mafft ${FILE} > ${FILE}.aligned
else
    echo "exception: Alignment method name seems wrong..."
fi

# Tree construction
# ML (RAxML)
if [ $TREEMETHOD = "raxmlML" ]; then
    time ${SOFTWARE} -s ${FILE}.aligned -n raxml -T $THREAD_NUM ${OPTION} -m ${MODEL} -p 12345
    wait
    mv RAxML_bestTree.raxml ${FILE}.aligned.tree
# NJ
elif [ $TREEMETHOD = "rapidnjNJ" ]; then
    time $SOFTWARE -c ${THREADNUM} ${FILE}.aligned -i fa ${OPTION} > ${FILE}.aligned.tree
# MP
elif [ $TREEMETHOD = "raxmlMP" ]; then
    time ${SOFTWARE} -s ${FILE}.aligned -y -n raxml -T $THREAD_NUM -m ${MODEL} -p 12345 # -y : only compute maximum parsimony!
    mv RAxML_parsimonyTree.raxml ${FILE}.aligned.tree
# ML (fasttree)
elif [ $TREEMETHOD = "fasttreeML" ]; then
    export OMP_NUM_THREADS=$THREADNUM
    time $SOFTWARE -gtr -nt ${OPTION} -seed 111 < ${FILE}.aligned > ${FILE}.aligned.tree
# others
else [ $TREEMETHOD = "unspecified" ]; then
    time $SOFTWARE ${FILE}.aligned > ${FILE}.aligned.tree
else
    echo "exception: Tree construction method name seems wrong..."
fi