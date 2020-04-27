echo '########################################################'
echo '##            Sample codes of FRACTAL 1.4             ##'
echo '##                                                    ##'
echo '##  Last update: April 24, 2020                       ##'
echo '##  Contact: Naoki Konno (naoki@bs.s.u-tokyo.ac.jp)   ##'
echo '##           Nozomu Yachie (nzmyachie@gmail.com)      ##'
echo '########################################################'
echo ''
echo '#####################################################################################################'
echo '# Please note that the deep distributed computing mode is disabled here because UGE is not          #'
echo '# supported by Code Ocean. The original code is avaialble at https://github.com/yachielab/FRACTAL.  #'
echo '#####################################################################################################'

chmod u+x /root/capsule/code/FRACTAL/FRACTAL.sh
export PATH=/root/capsule/code/FRACTAL:$PATH

echo -e "\nRunning /code/run.sh..."

echo "The expected total calculation time is 8-10 minutes."

# Example code 1
echo -e "\nRunning Example code 1..."
FRACTAL.sh -i FRACTAL/example/test.fa > FRACTALout.log
mv FRACTALout.nwk /results
mv FRACTALout.log /results
echo "FRACTALout.nwk and FRACTALout.log were generated in /results"

# Example code 2
echo -e "\nRunning Example code 2..."
FRACTAL.sh -i FRACTAL/example/test.fa -f FRACTAL_NJ -m rapidnjNJ > FRACTAL_NJ.log
mv FRACTAL_NJ.nwk /results
mv FRACTAL_NJ.log /results
echo "FRACTAL_NJ.nwk and FRACTAL_NJ.log were generated in /results"

# Example code 3
echo -e "\nRunning Example code 3..."
FRACTAL.sh -i FRACTAL/example/test.fa -f FRACTAL_ML -m fasttreeML -p "-fastest -quiet" -c 16 > FRACTAL_ML.log
mv FRACTAL_ML.nwk /results
mv FRACTAL_ML.log /results
echo "FRACTAL_ML.nwk and FRACTAL_ML.log were generated in /results"

echo -e "\nEnd of operation."