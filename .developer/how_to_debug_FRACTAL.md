### How to debug FRACTAL

0. I recommend to use Visual Studio Code with "Remote - SSH" extention as a developing platform

1. Activate the following line in FRACluster.new.py by the removing commentout mark `#`

    ```
    shutil.copytree(WD + "/INPUT", WD + "/INPUT_copy")
    ```

2. Execute with '-e' option: FRACTAL quits deleting intermediate files

3. See `FRACTALout/err` directory and identify the directory `dXX` where an error occured

4. Create a copy of the directory and initialize the directory by

    ```
    cp -r FRACTALout/nodes/dXX FRACTALout/nodes/dXX_copy 
    rm -r FRACTALout/nodes/dXX
    cp -r FRACTALout/nodes/dXX_copy/INPUT_copy FRACTALout/nodes/dXX/
    ```

5. Debug the code and retry the calculation for the directory by executing the corresponding script saved in `FRACTALout/executed`