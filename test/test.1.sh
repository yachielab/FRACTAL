for GZIP in "" "-g"; do
    for MODE in "-E" "-p ML" "-p MP" "-u -p ML" "-u -p MP"; do
        for DISTRIBUTION in "-d 1" "-d 2" "-d 100"; do
            string=$(echo $GZIP $MODE $DISTRIBUTION | tr -d ' ' | tr '-' )
            FRACTAL -i test.fa.gz $GZIP $MODE $DISTRIBUTION -j $string -f $string &
        done
    done 
done