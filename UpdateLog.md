### Update Log

<h4> 05/05/2019 Version 0.1.6 </h4>

1. a lot of bug fix
2. first version for benchmark

<h4> 02/25/2019 Version 0.1.5 </h4>

1. code update: enable distribution of phylogenetic placement.
2. bug fix: use user's own path even in calculation node

<h4> 02/04/2019 Version 0.1.4 </h4>
<h4> important! This version is following Version 0.1.2, (not 0.1.3) </h4>

1. code update: Implement options of FRACTAL.sh / Tree.sh <br>
2. code update: Enabled FastTreeMP

<h4> 02/04/2019 Version 0.1.3 </h4>

1. code update: Implement rooting method by phylogenetic placement.

<h4> 01/25/2019 Version 0.1.2 </h4>

1. code update : enabled ML (raxml), NJ (rapidnj) as well as MP (parsimonator in raxml) <br>
2. code update : tested in reedbush <br>
3. bug fixed : iterative sampling <br>

<h4> 01/20/2019 Version 0.1.1 </h4>

1. enabled qsub of FRACTAL-0.1.1

<h4> 12/24/2018 Version 0.1.0 </h4>

1. implementation : quit using ancestral reconstruction for assignment, but use ML-based method EPA-ng. <br>
2. implementation : quit choosing 1 subsample after multiple subsampling, but do itrative subsampling. <br>
3. implementation : offer multiple rooting method <br>
4. NOTICE : no qsub occurs in this version, because DDBJ supercomputer is busy. <br>
5. NOTICE : even though no qsub occurs, this version do multi-threading when EPA-ng works. <br>
         

<h4> 10/23/2018 Version 0.0.1 </h4>

1. implementation : mapping onto subsample tree = “MapOnTree”, and enabled the choice of "Mix"  <br>
2. code update : enable users to set number of computers & implemented no-qsub version <br>
3. code update : measure p-distance between IUPAC ambiguity sequences <br>

<h4> 10/17/2018 Version 0.0.0 </h4>

1. implementation : multiple subsampling was implemented. <br>
2. implementation : interactive user interface <br>
3. code update : halting condition of "qsubing" <br>