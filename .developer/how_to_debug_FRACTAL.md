### How to debug FRACTAL

0. I recommend to use Visual Studio Code with "Remote - SSH" extention as a developing platform

1. Activate the following line in FRACluster.new.py by the removing commentout mark `#`

    ```
    shutil.copytree(WD + "/INPUT", WD + "/INPUT_copy")
    ```

2. Execute FRACTAL with '-e' option: FRACTAL quits deleting intermediate files

3. See `FRACTALout/err` directory and identify the directory `FRACTALout/nodes/dXX` which caused the error 

4. Create a copy of the directory and initialize the directory by

    ```
    cp -r FRACTALout/nodes/dXX FRACTALout/nodes/dXX_copy 
    rm -r FRACTALout/nodes/dXX
    cp -r FRACTALout/nodes/dXX_copy/INPUT_copy FRACTALout/nodes/dXX/
    ```

5. Debug the code and retry the calculation only for the directory by executing the corresponding script saved in `FRACTALout/executed`

6. Execute a test script to check if FRACTAL can work for any combination of options 

    ```
    cd FRACTALout/test
    bash test.1.sh
    ```