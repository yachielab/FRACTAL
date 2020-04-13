FASTA_file_path=$1
raxmlHPC-SSE3 -s ${FASTA_file_path} -n test -m GTRGAMMA -p 1 > /dev/null
mv test_bestTree.raxml ${FASTA_file_path}.tree