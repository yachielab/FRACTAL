<h2>FRACTAL Installation and User Manual</h2>


### Overview of FRACTAL


**FRACTAL** (framework for distributed computing to trace large accurate lineages) is a new deep distributed computing framework that is designated to reconstruct extremely large lineages of nucleotide sequences using a software tool of choice. In brief, FRACTAL first subsamples a small number of sequences to reconstruct only an upper hierarchy of a target lineage and assigns the remaining sequences to its downstream clades, for each of which the same procedure is recursively iterated. Since the iteration procedure can be performed in parallel in a distributed computing framework, FRACTAL allows highly scalable reconstruction of an extremely large lineage tree.

<img src=images/fractal_concept.jpg width=10000x3000>


**Figure 1. Schematic diagram of FRACTAL.** From an input sequence pool, a given number of sequences are first randomly subsampled (Step 1) and their sample lineage tree is reconstructed with a rooting or provisional rooting sequence by lineage estimation software of choice (Step 2). Note that while an outgroup is preferred as the rooting sequence if available, any sequence can be used for unrooted tree estimation. The remaining input sequences are then mapped to the branches of the sample tree by phylogenetic placement (Step 3). If all of the input sequences are mapped on downstream branches of the sample tree so as to separate them into multiple distinct clades, their upstream lineage is considered to be true and fixed (Step 4), and the sequence group in each downstream clade is recursively subjected to the first process in a distributed computing node (Step 5). If any sequence(s) mapped on the root branch does not allow the grouping of input sequences into clades, the phylogenetic placement is repeated against a new sample tree generated for sequences randomly chosen from a union of the previous subsampled sequences and the “problematic” sequences (Step 6). This process generates a new sample tree in a biased manner such that it harbors the previous problematic sequences in its leaves and decreases the likelihood of acquiring problematic sequences in the following phylogenetic placement step. This procedure is repeated until the problem is solved, but only up to a given threshold number of times as long as the number of problematic sequences continues to be reduced in every retrial step. When the retrial cycle stops without solving the problem, the remaining problematic sequences are discarded, and the other sequence sets are separated into distinct clades and subjected to the first process. Accordingly, FRACTAL hierarchically generates expanding parallel computing trajectories, where each distributed computing job recursively generates a large set of successive jobs. When the number of input sequences is reduced to a certain threshold (hereafter called the naïve computing threshold) while the FRACTAL iteration cycles, the remaining marginal lineage is directly reconstructed by using the software of choice and the operation terminates for this computing trajectory (Step 7). For unrooted lineage estimation, the provisional rooting sequence can be removed after completion of the whole computation. Accordingly, FRACTAL enables efficient reconstruction of a large lineage by distributed computing while utilizing limited computing power and memory per node. FRACTAL is also effective even for a single computing node because its memory consumption level can be kept down for large lineage reconstructions.

### Supported Environment

1. FRACTAL can be executed on Linux OS
2. The distributed computing mode of FRACTAL requires UGE (Univa Grid Engine)

### Software Dependency

<h4>Required</h4>

1. Python3 (newer than 3.7.0) with Biopython module (required)
2. RAxML (raxmlHPC-PTHREADS-SSE3 and raxmlHPC-SSE3) (required)
3. EPA-ng (required)
4. RapidNJ (optional; if you want to use NJ for lineage reconstruction)
5. FastTreeMP (optional; if you want to use ML for lineage reconstruction)

### Software installation


#### Installation of FRACTAL


1. Download FRACTAL by

   ```shell
    git clone https://github.com/yachielab/FRACTAL.git
   ```

   or you can also obtain FRACTAL as follows

   ```shell
    wget https://github.com/yachielab/FRACTAL/archive/master.zip
    unzip master.zip
   ```

2. Add the absolute path of FRACTAL directory to $PATH
3. Make FRACTAL.sh executable

       chmod u+x FRACTAL.sh

#### Installation of [Anaconda](https://www.anaconda.com/distribution/) (required)

1. Execute the following commands

    ```
    wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
    bash Anaconda3-2018.12-Linux-x86_64.sh
    bash
    ```

2. Please set $PATH to anaconda3/bin.

#### Installation of [Biopython](https://anaconda.org/anaconda/biopython) (required)

1. Install Biopython by

   ```
   conda install biopython
   ```

#### Installation of [RAxML](https://github.com/stamatak/standard-RAxML) (required)

1. Execute the following commands

   ```
   wget https://github.com/stamatak/standard-RAxML/archive/master.zip
   unzip master.zip
   cd standard-RAxML-master
   make -f Makefile.SSE3.gcc
   rm *.o
   make -f Makefile.SSE3.PTHREADS.gcc
   rm *.o
   ```

2. Set $PATH to raxmlHPC-SSE3 & raxmlHPC-PTHREADS-SSE3 executable.

#### Installation of [EPA-ng](https://github.com/Pbdas/epa-ng) (required)

1. Install EPA-ng by

   ```
   conda install -c bioconda epa-ng
   ```

#### Installation of [RapidNJ](http://birc.au.dk/software/rapidnj/) (optional)

1. Install RapidNJ

   ```
   conda install -c bioconda rapidnj
   ```

#### Installation of [FastTreeMP](http://www.microbesonline.org/fasttree/) (optional)

1. Execute the following commands

   ```
   wget http://www.microbesonline.org/fasttree/FastTreeMP
   chmod u+x FastTreeMP
   ```

2. Please set $PATH to FastTreeMP executable.

### Sample Codes

The FRACTAL package contains an example input file in the `examples` directory so users can check the software functions as follows:

**Example 1**

Lineage estimation of the sequences in the test input file `test.fa` with the default parameter set without distributed computing.

```shell
FRACTAL.sh -i test.fa
```


Input

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa)  (FASTA format)

Output

​	 [`FRACTALout.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTALout.nwk) (Newick format) will be created in your current working directory

**Example 2**

Lineage estimation by NJ using RapidNJ with distributed computing where the maximum job number is set to 100. The output file name is set to `FRACTAL_NJ`.

```shell
FRACTAL.sh -i test.fa -f FRACTAL_NJ -m rapidnjNJ -d 100
```

Input:

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:

​	 [`FRACTAL_NJ.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_NJ.nwk) (Newick format) for intermediate files will be created in your working directory.

**Example 3**

Lineage estimation by ML using FastTreeMP with its option `-fastest -quiet` without distributed computing. The number of threads required for the phylogenetic placement and the sample tree reconstruction procedures is set to be 16. The output file name is set to `FRACTAL_ML`.  

```
FRACTAL.sh -i test.fa -f FRACTAL_ML -m fasttreeML -p "-fastest -quiet" -c 16
```

Input:    

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:  

​	 [`FRACTAL_ML.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_ML.nwk) (Newick format) 

### FRACTAL Usage

```
Usage:
    FRACTAL.sh
    [-v] [-h] [-i input_file] [-o output_file_path] [-f output_file_name]
    [-m method] [-p "options"] [-k sequence_number] [-b model_name]
    [-x iteration_number] [-t sequence_number]
    [-d job_number] [-c thread_number] [-e]
    [-r integer] [-O qsub_option] [-I first_qsub_option] [-j job_name]

Options:
    -v
      Print FRACTAL version; ignore all the other parameters
    -h
      Print the usage of FRACTAL; ignore all the other parameters
    -i <String>
      Input FASTA file
    -o <String>
      Output directory path. Default: current working directory
    -f <String>
      Output file name. Default: FRACTALout
    -s <String>
      Absolute path of the executable file used to reconstruct lineage tree in each iteration cycle.
    -p "<String>"
      Options for the software corresponding to the method selected by -m
    -k <Integer>
      Number of sequences for the subsampling procedure. Default: 100
    -b <String>
      Substitution model of RAxML for phylogenetic placement. Default: GTRCAT
    -x <Integer>
      Threshold for the maximum number of retrial iterations in the subsampling process
    -t <Integer>
      Threshold number of input sequences to switch to direct lineage tree reconstruction 
        in each iteration cycle. Default: 500
    -d <Integer>
      Maximum number of jobs permissible for distributed computing.
        Default: 1 (no distributed computing)
    -c <Integer>
      Number of threads for the subsample tree reconstruction and the phylogenetic placement
        in each iteration cycle. Default: 1
    -e
      Output intermediate files
    -r <Integer>
      Seed number for generation of random values. Default: 0
    -O "<String>"
      Options for qsub. Default: ""
      example:  -O "-pe def_slot 4 -l s_vmem=16G -l mem_req=16G" 
    -I "<String>"
      Options especially for the first qsub. Default: the string specified by -O
    -j "<String>"
      Name of the jobs distributed by FRACTAL. Default: "FRACTAL"
    -m <String, Permissible values: ‘raxmlMP’, ‘rapidnjNJ’ and ‘fasttreeML’>
      Method to reconstruct lineage tree in each iteration cycle. Default: raxmlMP
      When you specify -s option, this option will be ignored.
```


### Contact


Naoki Konno (The University of Tokyo) [naoki@bs.s.u-tokyo.ac.jp](mailto:naoki@bs.s.u-tokyo.ac.jp)

Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)