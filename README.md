<h2>FRACTAL Installation and User Manual</h2>


### Overview of FRACTAL


**FRACTAL** (framework for distributed computing to trace huge accurate lineage) is a new distributed computing framework that is designated to reconstruct a large lineage of DNA sequences that have evolved or diversified from an ancestral DNA sequence through the accumulation of mutations. FRACTAL enables an iterative top-down clustering of a large number of input sequences for lineage reconstruction by first randomly subsampling a defined number of sequences  in order to reconstruct a sample lineage tree.
The software then maps each initial input sequence on branches of the sample tree by phylogenetic placement. This result is used to group the input sequences into distinct clades with the information of the fixable upstream lineage tree for these clades. The same procedure is then recursively repeated to progressively reconstruct trees that are further downstream. Since the iteration procedure can be performed in parallel in a distributed computing framework, FRACTAL allows highly scalable reconstruction of an extremely large lineage tree.

<img src=images/fractal_concept.jpg width=10000x3000>


#### Figure. 1


Figure 1. Schematic diagram of FRACTAL. The software recursively iterates top-down clustering of input sequences for a large lineage tree reconstruction. (1) *s* number of sequences are randomly subsampled from a large input sequence pool. (2) The subsampled sequences are used to reconstruct a sample tree by any of the lineage reconstruction methods (the current version employs NJ: neighbor joining, MP: maximum parsimony, and ML: maximum likelihood). (3) Each input sequence is then mapped to branches of the sample tree by phylogenetic placement (the current version employs only ML). (4) If all of the input sequences are mapped on downstream branches of the sample tree, (5) they are separated into distinct clades, each of which is subjected to the same recursive procedure that can be performed in a distributed computing node. (6) If any sequence(s) is mapped on the root branch of the sample tree in (4) which does not allow the clustering of sequences into subgroups, the subsample procedure is repeated from a union of the previous subsampled sequences and the “problematic” sequence(s). This retrial procedure facilitates the generation of a new sample tree in a biased manner such that it tends to harbor the problematic sequences in its downstream branch(es). This retrial procedure (6) repeats up to *x* times as long as the reduction in number of problematic sequences is observed in every iterative step. When the retrial iteration stops, the remaining problematic sequences are discarded, and the other sequences mapped on the tree are subjected to the process (5). Finally, when the number of input sequences has become reduced in size to a certain point in the recursive processes running in parallel, the operation ends with the reconstruction of the marginal lineage for the remaining input sequences.

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
    git clone https://github.com/yachielab/FRACTAL
   ```

   or you can also obtain FRACTAL as follows

   ```shell
    wget ~~~~/FRACTAL_latest.zip
    unzip FRACTAL_latest.zip
   ```

2. Add the absolute path of FRACTAL directory to $PATH
3. Make FRACTAL.sh executable

       cd FRACTAL
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

​	 [`test.fa`]()  (FASTA format)

Output

​	 [`FRACTALout.nwk`]() (Newick format) will be created in your current working directory

**Example 2**

Lineage estimation by NJ using RapidNJ with distributed computing where the maximum job number is set to 100. The output file name is set to `FRACTAL_NJ`. Deletion of intermediate files is inactivated.

```shell
FRACTAL.sh -i test.fa -o FRACTAL_NJ -m rapidnjNJ -q 100
```

Input:

​	 [`test.fa`]() (FASTA format)

Output:

​	 [`FRACTAL_NJ.nwk`]() (Newick format) for intermediate files will be created in your working directory.

**Example 3**

Lineage estimation by ML using FastTreeMP with its option `-fastest -quiet` without distributed computing. The number of threads required for the phylogenetic placement and the sample tree reconstruction procedures is set to be 16. The output file name is set to `FRACTAL_ML`.  

```
FRACTAL.sh -i test.fa -o FRACTAL_ML -m fasttreeML -p "-fastest -quiet" -c 16
```

Input:    

​	 [`test.fa`]() (FASTA format)

Output:  

​	 [`FRACTAL_ML.nwk`]() (Newick format) 

### FRACTAL Usage

```
Usage:
    FRACTAL.sh
    [-v] [-h] [-i input_file] [-d output_file_path] [-o output_file_name]
    [-m method] [-p "options"] [-s sequence_number] [-b model_name]
    [-x iteration_number] [-l sequence_number]
    [-q job_number] [-t thread_number] [-e]
    [-r integer] [-M memory_size] [-I memory_size]

Options:
    -v
      Print FRACTAL version; ignore all the other parameters
    -h
      Print the usage of FRACTAL; ignore all the other parameters
    -i <String>
      Input FASTA file
    -d <String>
      Output directory path. Default: current working directory
    -o <String>
      Output file name. Default: FRACTALout
    -m <String, Permissible values: ‘raxmlMP’, ‘rapidnjNJ’ and ‘fasttreeML’>               
      Method to reconstruct lineage tree in each iteration cycle. Default: raxmlMP
    -p <”String”>
      Options for the software corresponding to the method selected by -m
    -s <Integer>
      Number of sequences for the subsampling procedure. Default: 100
    -b <String>
      Substitution model of RAxML for phylogenetic placement. Default: GTRCAT
    -x <Integer>
      Threshold for the maximum number of retrial iterations in the subsampling process
    -t <Integer>
      Threshold number of input sequences to switch to direct lineage tree reconstruction 
        in each iteration cycle. Default: 500
    -q <Integer>
      Maximum number of jobs permissible for distributed computing.
        Default: 1 (no distributed computing)
    -c <Integer>
      Number of threads for the subsample tree reconstruction and the phylogenetic placement
        in each iteration cycle. Default: 1
    -e
      Output intermediate files.
    -r <Integer>
      Seed number for generation of random values. Default: 0
    -M <Integer>GB
      Memory requirement per distributed computing node. Default: 16GB
    -I <Integer>GB
      Memory requirement for the first job in the distributed computing mode.
        Default: 16GB
```


### Note


For those of you who use the SHIROKANE supercomputer at the Human Genome Center of the University of Tokyo, please use the option "-k SHIROKANE". This will specify the OS to be used for FRACTAL to CentOS7.


### Contact


Naoki Konno (The University of Tokyo) [naoki@bs.s.u-tokyo.ac.jp](mailto:naoki@bs.s.u-tokyo.ac.jp)

Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)