<h2>FRACTAL Installation and User Manual</h2>

- [Overview of FRACTAL](#overview-of-fractal)
- [Supported Environment](#supported-environment)
- [Software Dependency](#software-dependency)
- [Software installation](#software-installation)
- [Sample Codes](#sample-codes)
- [FRACTAL Usage](#fractal-usage)

### Overview of FRACTAL

**FRACTAL** (framework for distributed computing to trace large accurate lineages) is a new deep distributed computing framework that is designated to reconstruct extremely large lineages of nucleotide sequences using a software tool of choice. In brief, FRACTAL first subsamples a small number of sequences to reconstruct only an upper hierarchy of a target lineage and assigns the remaining sequences to its downstream clades, for each of which the same procedure is recursively iterated (Fig. 1). Since the iteration procedures can be performed in parallel in a distributed computing framework, FRACTAL allows highly scalable reconstruction of an extremely large lineage tree. 
The input sequences of FRACTAL need not be unaligned (Fig. 2). When many unaligned sequences are queried to FRACTAL for lineage reconstruction, only the sequences subsampled from the input sequences are aligned by MSA for the sample tree reconstruction in each cycle. For phylogenetic placement, each of the remaining input sequences is mapped to the sample tree through “plus-one” alignment to the MSA result that reconstructed the sample tree. This strategy enables the scalable lineage reconstruction of unaligned sequences.

<img src=images/fractal_concept.jpg width=10000x3000>

**Figure 1. Schematic diagram of FRACTAL.** From an input sequence pool, a given number of sequences are first randomly subsampled (Step 1), and their sample lineage tree is reconstructed with a rooting or provisional rooting sequence by lineage estimation software of choice (Step 2). (Note that while an outgroup or a true ancestral sequence (for cell lineage tracing) is preferred as the rooting sequence if available, any sequence can be used for unrooted tree estimation.) Each of the remaining input sequences is then mapped to its most proximal branch of the sample tree by phylogenetic placement (Step 3). If all of the input sequences are mapped on downstream branches of the sample tree to separate them into multiple distinct clades, their upstream lineage is considered to be resolved (Step 4), and the sequence group in each downstream clade is recursively subjected to the first process in a distributed computing node (Step 5). If any sequence(s) is mapped on the root branch, which does not allow the grouping of input sequences into clades, the phylogenetic placement is repeated against a new sample tree generated for sequences randomly chosen from a union of the previous subsampled sequences and the “problematic” sequences (Step 6). This process generates a new sample tree in a biased manner such that it harbors the previous problematic sequences in its leaves and decreases the probability of acquiring problematic sequences in the following phylogenetic placement step. The retrial process is repeated until the problem is solved, but only up to a given threshold number of times and as long as the number of problematic sequences continues to be reduced in every retrial step. When the retrial cycle stops without completely solving the problem, the remaining problematic sequences are discarded, and the other sequence sets are separated into distinct clades and subjected to the first process in separate computing nodes. Accordingly, FRACTAL generates a hierarchy of expanding parallel computing trajectories, where each distributed computing job recursively generates a large set of successive jobs. When the number of input sequences is reduced to a certain threshold (hereafter called the naïve computing threshold), the remaining lineage is directly reconstructed using the designated software, and the operation terminates for this computing trajectory (Step 7). For unrooted lineage estimation, the provisional rooting sequence can be removed after completing the whole computation. FRACTAL, therefore, enables efficient reconstruction of large lineages by distributed computing while utilizing limited computing power and memory per node. FRACTAL is also effective even for a single computing node because its memory consumption level can be kept down even for large lineage reconstructions.

<img src=images/fractal_unaligned.jpg width=10000x3000>

**Figure 2. Workflow for unaligned datasets.** When FRACTAL is executed with unaligned sequences, their lineage is reconstructed with an incremental sequence alignment strategy as follows. For the first sample tree construction, sequences subsampled from the initial sequence pool are aligned by MSA (Step 2). For phylogenetic placement, an HMM profile is constructed from the MSA result, and each of the input sequences is independently aligned to the MSA of the subsampled sequences (Step 4i). For each plus-one alignment result, after removing nucleotides that are not matched to the reference MSA result, the input sequence is subjected to phylogenetic placement to determine its most proximal branch of the sample tree (Step 4ii). To save the computing cost, the MSA and plus-one alignment results are inherited to the next FRACTAL cycles certain times, in which the MSA and plus-one alignment results of sequences are reused according to the alignment coordinates kept from the parental results. Assuming overall similarity of input sequences in each FRACTAL cycle elevates along with the top-down lineage reconstruction process, the MSA for sample tree reconstruction and plus-one alignment results are updated when the input sequence pool size of a FRACTAL iteration is reduced to a fold-change threshold or less relative to that of the last update. When the size of input sequence pool is less than the naïve computing threshold, MSA is performed for the entire input sequences followed by terminal lineage reconstruction. 

### Supported Environment

1. FRACTAL can be executed on Linux OS
2. The distributed computing mode of FRACTAL requires UGE (Univa Grid Engine)

### Software Dependency

1. Python3 (version: 3.7.0 or later) with Biopython (version: 1.76) module *required*
2. RAxML (raxmlHPC-PTHREADS-SSE3 and raxmlHPC-SSE3) (version: 8.2.12) *required*
3. EPA-ng (version: 0.3.5) *required*
4. Seqkit (version: 0.11.0) *required*
5. Trimal (version: 1.4.1) *required*
6. RapidNJ (version: 2.3.2) *optional; if you want to use NJ for lineage reconstruction*
7. FastTreeMP (version: 2.1.10) *optional; if you want to use ML for lineage reconstruction*
8. MAFFT (version: 7.464) *optional; if you want to reconstruct a lineage from unaligned sequences*
9. HMMER (version: 3.3) *optional; if you want to reconstruct a lineage from unaligned sequences*


### Software installation

Each installation step will take less than ~1 min


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

2. Please set `$PATH` to `anaconda3/bin`.

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

### Sample Codes

The FRACTAL package contains an example input file in the `examples` directory so users can check the software functions as follows:

**Example 1**

Lineage estimation of the sequences in the test input file `test.fa` with the default parameter set without distributed computing. The computation will take several minutes.

```shell
FRACTAL -i test.fa
```


Input:

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa)  (FASTA format)

Output:

​	 [`FRACTALout.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTALout.nwk) (Newick format) will be created in your current working directory

**Example 2**

Lineage estimation by sample tree estimation using RapidNJ with distributed computing where the maximum number of computing nodes is set to 10. The output file name is set to `FRACTAL_NJ`. The computation will take several minutes.

```shell
FRACTAL -i test.fa -f FRACTAL_NJ -m rapidnjNJ -d 10
```

Input:

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:

​	 [`FRACTAL_NJ.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_NJ.nwk) (Newick format) for intermediate files will be created in your working directory.

**Example 3**

Lineage estimation by sample tree estimation using FastTreeMP with its option `-fastest -quiet` without distributed computing. The number of threads required for the phylogenetic placement and the sample tree reconstruction procedures is set to be 16. The output file name is set to `FRACTAL_ML`.  The computation will take ~5 min.

```shell
FRACTAL -i test.fa -f FRACTAL_ML -m fasttreeML -a "-fastest -quiet" -c 16
```

Input:    

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:  

​	 [`FRACTAL_ML.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_ML.nwk) (Newick format) 

**Example 4**

Lineage estimation with a software tool of choice and user defined parameters.

1. Prepare a shell script that takes a FASTA file as an input file, calculate a lineage of the input sequences, and output it to a Newick format file whose name is inherited from the input FASTA file, like from `foo.fa` to `foo.fa.tree`.

2. Add the absolute path of the shell script to `$PATH`

3. Make the shell script executable

4. Execute FRACTAL as follows. The example shell script file [`ml_raxml.sh`](https://github.com/yachielab/FRACTAL/blob/hotfix/example/script/ml_raxml.sh) below is prepared and provided in the installation package for ML method with GTR-Gamma model by RAxML. The maximum number of computing nodes is set to 10 in the following command. The computation will take ~5 min.

   ```shell
   FRACTAL -i test.fa -f FRACTAL_raxml -s ml_raxml.sh -d 10
   ```

Input:    

​	 [`test.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.fa) (FASTA format)

Output:  

​	 [`FRACTAL_raxml.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_raxml.nwk) (Newick format) 

**Example 5**

Lineage estimation from an unaligned sequence dataset. The number of sequences in a subsample and the threshold number of input sequences to swich to direct lineage computing are set to be 1,000.

```shell
FRACTAL -i test.unaligned.fa -f FRACTAL_unaligned -k 10 -t 10 -m fasttreeML -u
```

Input:    

​	 [`test.unaligned.fa`](https://github.com/yachielab/FRACTAL/blob/master/example/test.unaligend.fa) (FASTA format)

Output:  

​	 [`FRACTAL_unaligned.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_unaligned.nwk) (Newick format) 

**Example 6**

Lineage estimation from sets of substitutions, insertions and deletions  using MP (RAxML) and phylogenetic placement using MP (RAxML). When `-E` option is specified, FRACTAL takes a special input file format to describe any set of mutations (see example input).

```shell
FRACTAL -i test.edit -f FRACTAL_edit -p MP -E
```

Input:    

​	 [`test.edit`](https://github.com/yachielab/FRACTAL/blob/master/example/test.edit) (Original format: One sequence is represented in one row as "\<Name of a sequence\>\t\<Mutation pattern name 1\>;\<Mutation pattern name 2\>;...")

Output:  

​	 [`FRACTAL_edit.nwk`](https://github.com/yachielab/FRACTAL/blob/master/example/output/FRACTAL_edit.nwk) (Newick format) 

> Note for option -E:
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
      Print FRACTAL version; ignore all the other parameters
    -h
      Print the usage of FRACTAL; ignore all the other parameters
    -i <String>
      Input FASTA file
    -o <String>
      Output directory path. Default: current working directory
    -f <String>
      Output file name. Default: FRACTALout
    -m <String, Permissible values: ‘raxmlMP’, ‘rapidnjNJ’ and ‘fasttreeML’>
      Method to reconstruct lineage tree in each iteration cycle. Default: raxmlMP
        When you specify -s option, this option will be ignored.
    -a "<String>"
      Options for the software corresponding to the method selected by -m
    -s <String>
      File name of a shell script used to reconstruct lineage tree in each iteration cycle.
        See sample codes (example 4).
    -k <Integer>
      Number of sequences for the subsampling procedure. Default: 10
    -u
      Reconstuct a lineage from unaligned sequences.
    -w <float>
      Fold-change threshold which determines the timing to update the MSA and the plus-one 
        alignment result. This option is valid only when -u is specified. Default: 0.5
    -b <String>
      Substitution model of RAxML for phylogenetic placement. Default: GTRCAT
    -p <String, Permissible values: ‘ML’, ‘MP’>
      Method of phylogenetic placement Default: ML
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
      Reconstruct a lineage from lists of mutation patterns. 
      Currently, this option can be specified only with "-p MP -d 1", and without "-g".
    -O "<String>"
      Options for qsub. Default: ""
        example:  -O "-pe def_slot 4 -l s_vmem=16G -l mem_req=16G" 
    -I "<String>"
      Options especially for the first qsub. Default: the string specified by -O
    -A "<String>"
      Options especially for the last qsub (tree assembly). 
        Default: the string specified by -O
    -j "<String>"
      Name of the jobs distributed by FRACTAL. Default: "FRACTAL"
    -l <Integer>
      Maximum number of FRACTAL iterations. Default: 10000
    -z <Integer>
      Number of tips extracted from the sample tree. 
        Default: the same value as the value specified by -k
    -g 
      Gzip intermediate files.
```


### Contact


Naoki Konno (The University of Tokyo) [konno-naoki555@g.ecc.u-tokyo.ac.jp](mailto:konno-naoki555@g.ecc.u-tokyo.ac.jp)

Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)

