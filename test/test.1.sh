WD=$(pwd)

if [ -e test ]; then
    rm -r test
fi
mkdir test

idx=0

# for aligned FASTA file 
for GZIP in "" "-g"; do
    for MODE in "-p ML" "-p MP"; do
        for DISTRIBUTION in "-d 1" "-d 2" "-d 100"; do
            string=a$(echo $GZIP $MODE $DISTRIBUTION | tr -d ' ' | tr -d '-' )
            mkdir $WD/test/$string; cd $WD/test/$string
            (FRACTAL -i ${WD}/test.fa $GZIP $MODE $DISTRIBUTION -j _${idx}_ -f test${idx}_$string -e &> $string.outerr 
            echo "$WD/test/$string/test${idx}_${string}.nwk"               >> $WD/test/$string/result.out
            tree_info --Ntips -t $WD/test/$string/test${idx}_${string}.nwk >> $WD/test/$string/result.out) &

            idx=$(expr $idx + 1)
        done
    done 
done

# for unaligned FASTA file 
for GZIP in "" "-g"; do
    for MODE in "-u -p ML" "-u -p MP"; do
        for DISTRIBUTION in "-d 1" "-d 2" "-d 100"; do
            string=u$(echo $GZIP $MODE $DISTRIBUTION | tr -d ' ' | tr -d '-' )
            mkdir $WD/test/$string; cd $WD/test/$string
            (FRACTAL -i ${WD}/test.unaligned.fa.gz $GZIP $MODE $DISTRIBUTION -j _${idx}_ -f test${idx}_$string -e &> $string.outerr 
            echo "$WD/test/$string/test${idx}_${string}.nwk"               >> $WD/test/$string/result.out
            tree_info --Ntips -t $WD/test/$string/test${idx}_${string}.nwk >> $WD/test/$string/result.out) &

            idx=$(expr $idx + 1)
        done
    done 
done


# for edit list file 
for GZIP in ""; do
    for MODE in "-E -p MP"; do
        for DISTRIBUTION in "-d 1"; do
            string=i$(echo $GZIP $MODE $DISTRIBUTION | tr -d ' ' | tr -d '-' )
            mkdir $WD/test/$string; cd $WD/test/$string
            (FRACTAL -i ${WD}/test.indel.gz $GZIP $MODE $DISTRIBUTION -j _${idx}_ -f test${idx}_$string -e -k 20 -t 100 &> $string.outerr 
            echo "$WD/test/$string/test${idx}_${string}.nwk"               >> $WD/test/$string/result.out
            tree_info --Ntips -t $WD/test/$string/test${idx}_${string}.nwk >> $WD/test/$string/result.out) &

            idx=$(expr $idx + 1)
        done
    done 
done

# to use custom tree reconstruction software
for GZIP in "-g"; do
    for MODE in "-p ML"; do
        for DISTRIBUTION in "-d 100"; do
            string=s$(echo $GZIP $MODE $DISTRIBUTION | tr -d ' ' | tr -d '-' )
            mkdir $WD/test/$string; cd $WD/test/$string
            (FRACTAL -i ${WD}/test.fa $GZIP $MODE $DISTRIBUTION -j _${idx}_ -f test${idx}_$string -s ml_raxml.sh -d 100 &> $string.outerr 
            echo "$WD/test/$string/test${idx}_${string}.nwk"               >> $WD/test/$string/result.out
            tree_info --Ntips -t $WD/test/$string/test${idx}_${string}.nwk >> $WD/test/$string/result.out) &

            idx=$(expr $idx + 1)
        done
    done 
done