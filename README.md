<h2>FRACTAL Installation and User Manual</h2>

- [Overview of FRACTAL](#overview-of-fractal)
- [Supported Environment](#supported-environment)
- [Software Dependency](#software-dependency)
- [Software installation](#software-installation)
- [Sample Codes](#sample-codes)
- [FRACTAL Usage](#fractal-usage)

### Overview of FRACTAL

**FRACTAL** (framework for distributed computing to trace large accurate lineages) is a deep distributed computing framework that reconstruct extremely large lineages of nucleotide sequences using a software tool of choice (Figure 1). From an input sequence pool, a given number of sequences are first randomly subsampled (Step 1). Their sample lineage tree is reconstructed with a rooting or provisional rooting sequence by lineage estimation software of choice (Step 2). (Note that while an outgroup or a true ancestral sequence, for cell lineage tracing, is preferred as the rooting sequence if available, any sequence can be used for unrooted tree estimation.) Each of the remaining input sequences is then mapped to its most proximal branch of the sample tree by phylogenetic placement (Step 3). If all of the input sequences are mapped on the downstream branches of the sample tree to separate them into multiple distinct clades, their upstream lineage is considered to be resolved (Step 4), and the sequence group in each downstream clade is recursively subjected to the first process in a distributed computing node (Step 5). If any sequences is mapped on the root branch, which does not allow the grouping of input sequences into clades, the phylogenetic placement is repeated against a new sample tree generated for sequences randomly chosen from a union of the previous subsampled sequences and the “problematic” sequences (Step 6). This process generates a new sample tree in a biased manner such that it harbors the previous problematic sequences in its leaves and decreases the probability of acquiring problematic sequences in the subsequent phylogenetic placement step. The retrial process is repeated until the problem is solved, but only up to a given threshold number of times and as long as the number of problematic sequences continues to be reduced in every retrial step. When the retrial cycle stops without completely solving the problem, the remaining problematic sequences are discarded, and the other sequence sets are separated into distinct clades and subjected to the first process in separate computing nodes. Accordingly, FRACTAL generates a hierarchy of expanding parallel computing trajectories, where each distributed computing job recursively generates a large set of successive jobs. When the number of input sequences is reduced to a certain threshold (hereafter called the naïve computing threshold), the remaining lineage is directly reconstructed using the designated software, and the operation terminates for this computing trajectory (Step 7). For unrooted lineage estimation, the provisional rooting sequence can be removed after completing the entire computation. Thus, FRACTAL enables the efficient reconstruction of large lineages by distributed computing while utilizing limited computing power and memory per node.

<img src=images/fractal_concept.jpg width=10000x3000>

**Figure 1. Schematic diagram of FRACTAL.** 

Multiple sequence alignment (MSA) of the input sequences is a prerequisite for some steps in FRACTAL (Figure 2). When unaligned sequences are queried to FRACTAL for lineage reconstruction, only the sequences subsampled from the input sequences are aligned by MSA for the sample tree reconstruction in each cycle. For phylogenetic placement, each of the remaining input sequences is mapped to the sample tree through “plus-one” alignment to the MSA result that reconstructed the sample tree. This strategy enables the scalable lineage reconstruction of unaligned sequences.

<img src=images/fractal_unaligned.jpg width=10000x3000>

**Figure 2. Evolutionary lineage reconstruction without preliminary MSA.** A given number of sequences are first randomly subsampled from the input sequences (Step 1). The subsampled sequences are aligned with a common root sequence by MSA using MAFFT (Step 2) and a sample tree is reconstructed by a software tool of choice (Step 3). Each of the remaining input sequences are then independently added to the MSA result by “plus-one” alignment using HMMER (Step 4i) and placed on the sample tree (Step 4ii).

### Supported Environment

1. FRACTAL can be executed on Linux OS
2. The distributed computing mode of FRACTAL requires UGE (Univa Grid Engine)

### Software Dependency

1. Python3 (version: 3.7.0 or later) with Biopython (version: 1.76) module *required*
2. RAxML (raxmlHPC-PTHREADS-SSE3 and raxmlHPC-SSE3) (version: 8.2.12) *required*
3. EPA-ng (version: 0.3.5) *required*
4. Seqkit (version: 0.11.0) *required*
5. Trimal (version: 1.4.1) *required*
6.	RapidNJ (version: 2.3.2) optional; for NJ method
7.	FastTreeMP (version: 2.1.10) optional; for ML method
8.	MAFFT (version: 7.464) optional; for lineage reconstruction with the incremental MSA method
9.	HMMER (version: 3.3) optional; for lineage reconstruction with the incremental MSA method


### Software installation

Each installation step will take less than 1 min


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

2. Add the absolute path of `FRACTAL` directory to `$PATH`

3. Make `FRACTAL` executable

   ```shell
   chmod u+x FRACTAL
   ```

#### Installation of [Anaconda](https://www.anaconda.com/distribution/) (required)

1. Execute the following commands

    ```shell
    wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
    bash Anaconda3-2018.12-Linux-x86_64.sh
    bash
    ```

2. Set `$PATH` to `anaconda3/bin`.

#### Installation of [Biopython](https://anaconda.org/anaconda/biopython) (required)

1. Install Biopython by

   ```shell
   conda install biopython
   ```

#### Installation of [RAxML](https://github.com/stamatak/standard-RAxML) (required)

1. Execute the following commands

   ```shell
   wget https://github.com/stamatak/standard-RAxML/archive/master.zip
   unzip master.zip
   cd standard-RAxML-master
   make -f Makefile.SSE3.gcc
   rm *.o
   make -f Makefile.SSE3.PTHREADS.gcc
   rm *.o
   ```

2. Set `$PATH` to `raxmlHPC-SSE3` & `raxmlHPC-PTHREADS-SSE3` executable.

#### Installation of [EPA-ng](https://github.com/Pbdas/epa-ng) (required)

1. Install EPA-ng by

   ```shell
   conda install -c bioconda epa-ng
   ```

#### Installation of [Seqkit](https://bioinf.shenwei.me/seqkit/) (required)

1. Install Seqkit by

   ```shell
   conda install -c bioconda seqkit
   ```

#### Installation of [RapidNJ](http://birc.au.dk/software/rapidnj/) (optional)

1. Install RapidNJ

   ```
   conda install -c bioconda rapidnj
   ```

#### Installation of [FastTreeMP](http://www.microbesonline.org/fasttree/) (optional)

1. Execute the following commands

   ```shell
   wget http://www.microbesonline.org/fasttree/FastTreeMP
   chmod u+x FastTreeMP
   ```

2. Please set `$PATH` to `FastTreeMP` executable.

#### Installation of [MAFFT](https://mafft.cbrc.jp/alignment/software/) (optional)

1. Install MAFFT

   ```shell
   conda install -c bioconda mafft
   ```

#### Installation of [HMMER](http://hmmer.org/) (optional)

1. Install HMMER

   ```shell
   conda install -c bioconda hmmer
   ```

#### Installation of [trimAl](http://trimal.cgenomics.org/) (optional)

1. Install trimAl

   ```shell
   conda install -c bioconda trimal
   ```

### Sample Commands

The FRACTAL package contains an example input file in the `examples` directory so users can check the software functions as follows:

**Example 1**

Lineage estimation of the sequences in the test input file `test.fa` with the default parameter set without distributed computing. This will take several minutes.

```shell
FRACTAL -i test.fa
```


Input:

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa)  (FASTA format)

Output:

​	 [`FRACTALout.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTALout.nwk) (Newick format) will be created in your current working directory

**Example 2**

Lineage estimation by sample tree estimation using RapidNJ with distributed computing where the maximum number of computing nodes is set to 10. The output file name is set to `FRACTAL_NJ`. This will take several minutes.

```shell
FRACTAL -i test.fa -f FRACTAL_NJ -m rapidnjNJ -d 10
```

Input:

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:

​	 [`FRACTAL_NJ.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_NJ.nwk) (Newick format) for intermediate files will be created in your working directory.

**Example 3**

Lineage estimation by sample tree estimation using FastTreeMP with its option `-fastest -quiet` without distributed computing. The number of threads required for the phylogenetic placement and the sample tree reconstruction procedures is set to be 16. The output file name is set to `FRACTAL_ML`. This will take approximately 5 minutes.

```shell
FRACTAL -i test.fa -f FRACTAL_ML -m fasttreeML -a "-fastest -quiet" -c 16
```

Input:    

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:  

​	 [`FRACTAL_ML.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_ML.nwk) (Newick format) 

**Example 4**

Lineage estimation with a software tool of choice and user defined parameters.

1.	Prepare a shell script that takes a FASTA format file and output their lineage to a Newick format file whose name is inherited from the input file name, like from foo.fa to foo.fa.tree.

2.	Add the absolute path of the shell script to $PATH

3.	Make the shell script executable

4. Execute FRACTAL as follows. The example shell script file [`ml_raxml.sh`](https://github.com/yachielab/FRACTAL/blob/hotfix/example/script/ml_raxml.sh) below is prepared and provided in the installation package for ML method with GTR-Gamma model by RAxML. The maximum number of computing nodes is set to 10 in the following command. This will take approximately 5 minutes.

   ```shell
   FRACTAL -i test.fa -f FRACTAL_raxml -s ml_raxml.sh -d 10
   ```

Input:    

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:  

​	 [`FRACTAL_raxml.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_raxml.nwk) (Newick format) 

**Example 5**

Lineage estimation without the preliminary MSA of the input sequences. The number of sequences for the subsampling process and the naïve computing threshold are both set to 10.

```shell
FRACTAL -i test.unaligned.fa -f FRACTAL_unaligned -k 10 -t 10 -m fasttreeML -u
```

Input:    

​	 [`test.unaligned.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.unaligned.fa) (FASTA format)

Output:  

​	 [`FRACTAL_unaligned.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_unaligned.nwk) (Newick format) 

**Example 6**

Lineage estimation from mutation lists of the input sequences using MP method for the phylogenetic placement.

```shell
FRACTAL -i test.edit -f FRACTAL_edit -p MP -E
```

Input:    

​	 [`test.edit`](https://github.com/yachielab/FRACTAL/blob/master/example/test.edit) (Each sequence is represented with a list of its mutations in a single row as follows)

```txt
[Sequence name]\t[Mutation pattern name 1];[Mutation pattern name 2];...
```

Output:  

​	 [`FRACTAL_edit.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_edit.nwk) (Newick format) 

**Note**:

While most of the practical phylogeny estimation tools can consider only substitutions, FRACTAL allows any software of choice to use both substitutions and indels for scalable lineage reconstruction. By specifying option -E, mutation patterns in the input sequences are converted into another sequence representation using four nucleotide sequences in the sample tree reconstruction and phylogenetic placement of each cycle as follows. All of the mutation patterns in sequences subsampled at each FRACTAL iteration are first encoded by binary letter patterns of unique positions to convert each mutated source sequence into a binary sequence of a fixed length, each of which represents the presence and absence of a specific mutation observed in the original subsampled pool by 1 and 0, respectively. The binary sequence of 1/0 is then further converted into a nucleotide sequence of thymine/cytosine or adenine/guanine at each position. The converted sequences are used to construct a sample tree using the given software of choice. Phylogenetic placement is achieved by converting the remaining sequences using the same encoding rule, ignoring the unique mutation patterns observed outside the original subsampled pool. The encoding rule is updated for every sample tree reconstruction, economically taking into account all mutations in the original input sequences through the lineage reconstruction process. Terminal tree reconstruction is performed using the same sequence conversion process. 


### FRACTAL Usage

```
Usage:
    FRACTAL.sh
    [-v] [-h] [-i input_file] [-o output_file_path] [-f output_file_name]
    [-m method] [-a "options"] [-s script_file_name] [-k sequence_number]
    [-u] [-w alignment_frequency] [-b model_name] 
    [-p placement_method] [-x iteration_number] [-t sequence_number]
    [-d job_number] [-c thread_number] [-e] [-r integer] [-n file_path] [-E] 
    [-O qsub_option] [-I first_qsub_option] [-A last_qsub_option] [-j job_name] 
    [-l iteration_upper_limit] [-g]

Options:
    -v
      Print FRACTAL version. Ignore all the other parameters.
    -h
      Print the usage of FRACTAL. Ignore all the other parameters.
    -i <String>
      Input FASTA file
    -o <String>
      Output directory path. Default: current working directory
    -f <String>
      Output file name. Default: FRACTALout
    -m <String, Permissible values: ‘raxmlMP’, ‘rapidnjNJ’ and ‘fasttreeML’>
      Method to reconstruct lineage tree in each iteration cycle. Default: raxmlMP
        When -s is specified, this option will be ignored.
    -a "<String>"
      Options for the software tool selected by -m
    -s <String>
      File name of a shell script for a user-defined lineage reconstruction tool 
    -k <Integer>
      Number of sequences for the subsampling procedure. Default: 10
    -z <Integer>
      Number of leaves randomly extracted from sample trees for the phylogenetic placement.
        Default: same as the value specified by -k
    -u
      Reconstuct a lineage from unaligned sequences.
    -w <float>
      Fold-change threshold that determines the timing to update the MSA and the plus-one 
        alignment result. This option is valid only when -u is specified. Default: 0.5
    -b <String>
      Substitution model of RAxML for phylogenetic placement. Default: GTRCAT
    -p <String, Permissible values: ‘ML’, ‘MP’>
      Method of phylogenetic placement. Default: ML
    -x <Integer>
      Threshold for the maximum number of retrial iterations in the subsampling process
    -t <Integer>
      Threshold number of input sequences to switch to direct lineage tree reconstruction 
        in each iteration cycle. Default: 10
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
    -E 
      Reconstruct a lineage from mutation lists of the input sequences. 
      Currently, this option can be specified only with "-p MP -d 1" and without "-g".
    -O "<String>"
      Options for qsub. Default: ""
        example:  -O "-pe def_slot 4 -l s_vmem=16G -l mem_req=16G" 
    -I "<String>"
      Options specified for the first qsub job. Default: Same as the string specified by -O
    -A "<String>"
      Options specified for the last qsub job (tree assembly). 
        Default: Same as the string specified by -O
    -j "<String>"
      Common name header of the jobs distributed by FRACTAL. Default: "FRACTAL"
    -l <Integer>
      Maximum number of FRACTAL iterations. Default: 10000
    -g 
      gzip intermediate files.
```


### Contact


Naoki Konno (The University of Tokyo) [konno-naoki555@g.ecc.u-tokyo.ac.jp](mailto:konno-naoki555@g.ecc.u-tokyo.ac.jp)

Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)

