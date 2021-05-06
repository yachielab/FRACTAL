echo '########################################################'
echo '##            Sample codes of FRACTAL 2.6             ##'
echo '##                                                    ##'
echo '##  Last update: April 03, 2021                       ##'
echo '##  Contact: Naoki Konno (naoki@bs.s.u-tokyo.ac.jp)   ##'
echo '##           Nozomu Yachie (nzmyachie@gmail.com)      ##'
echo '########################################################'
echo ''
echo '#####################################################################################################'
echo '# Please note that the deep distributed computing mode is disabled here because UGE is not          #'
echo '# supported by Code Ocean. The original code is avaialble at https://github.com/yachielab/FRACTAL.  #'
echo '#####################################################################################################'

chmod u+x /root/capsule/code/FRACTAL/FRACTAL
chmod u+x /root/capsule/code/FRACTAL/example/script/ml_raxml.sh
export PATH=/root/capsule/code/FRACTAL:/root/capsule/code/FRACTAL/example/script:$PATH

echo -e "\nRunning /code/run.sh..."

echo "The expected total calculation time is 15-20 minutes."

# Example code 1
echo -e "\nRunning Example code 1..."
FRACTAL -i FRACTAL/example/test.fa > FRACTALout.log
mv FRACTALout.nwk /results
mv FRACTALout.log /results
echo "FRACTALout.nwk and FRACTALout.log were generated in /results"


# Example code 2
echo -e "\nRunning Example code 2..."
FRACTAL -i FRACTAL/example/test.fa -f FRACTAL_NJ -m rapidnjNJ > FRACTAL_NJ.log
mv FRACTAL_NJ.nwk /results
mv FRACTAL_NJ.log /results
echo "FRACTAL_NJ.nwk and FRACTAL_NJ.log were generated in /results"


# Example code 3
echo -e "\nRunning Example code 3..."
FRACTAL -i FRACTAL/example/test.fa -f FRACTAL_ML -m fasttreeML -a "-fastest -quiet" -c 16 > FRACTAL_ML.log
mv FRACTAL_ML.nwk /results
mv FRACTAL_ML.log /results
echo "FRACTAL_ML.nwk and FRACTAL_ML.log were generated in /results"


# Example code 4
echo -e "\nRunning Example code 4..."
FRACTAL -i FRACTAL/example/test.fa -f FRACTAL_raxml -s /root/capsule/code/FRACTAL/example/script/ml_raxml.sh > FRACTAL_raxml.log
mv FRACTAL_raxml.nwk /results
mv FRACTAL_raxml.log /results
echo "FRACTAL_raxml.nwk and FRACTAL_raxml.log were generated in /results"


# Example code 5
echo -e "\nRunning Example code 5..."
FRACTAL -i FRACTAL/example/test.unaligned.fa -f FRACTAL_unaligned -k 10 -t 10 -m fasttreeML -u > FRACTAL_unaligned.log
mv FRACTAL_unaligned.nwk /results
mv FRACTAL_unaligned.log /results
echo "FRACTAL_unaligned.nwk and FRACTAL_unaligned.log were generated in /results"


# Example code 6
echo -e "\nRunning Example code 6..."
FRACTAL -i FRACTAL/example/test.edit -f FRACTAL_edit -p MP -E > FRACTAL_edit.log
mv FRACTAL_edit.nwk /results
mv FRACTAL_edit.log /results
echo "FRACTAL_edit_pattern.nwk and FRACTAL_edit_pattern.log were generated in /results"

echo -e "\nEnd of operation."