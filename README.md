# Pangenome graphs and their application to biodiversity genomics

This repository contains the script and files for the pratical part of the course "Pangenome graphs and their application to biodiversity genomics" held the 8th of September 2025 in Ferrara, Italy, as part of the [SIBE summer school](https://sites.google.com/view/sibesummerschool/home-page).


## Introduction

This course will...

## Table of content

0. First steps
    1. Files and folders
    2. Tools
1. Pangenome construction
    1. Fasta input files
    2. bTaeGut.seqfile input file
    
## 0. First steps

### 0.1 Files and folders

Copy this repository in your folder:

```
git clone https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics.git
```

The main directory of this repository contains all the input files needed for the exercises, as well as empty folders where your outputs will be written. All commands can be run from the main directory and they will output the data in the correct folders automatically.

A separate folder, [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data), contains all the expected outputs organized using the same folder structure. You can use these reference files to compare results with the expected outcome, complete exercises after the session at home and troubleshoot and debug each step.

### 0.2 Tools

Almost all the tools needed for this course are inside a conda environment.
Remember to **activate the conda environment** before running any command:

```
conda activate pangenomics
```
Other tools can be found here: 

## 1. Pangenome construction

To construct the pangenome we will use the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)<sup>1</sup>

1. A user-selected reference genome is used as the initial backbone
2. The reference is progressively augmented with structural variation from the other genomes by minigraph, a sequence-to-graph aligner, as a graph constructor. The resulting graph is SV only (>50 bp)
3. All assemblies are aligned back to the graph with a minimap2-like algorithm that generates base-level alignments for each reference chromosome separately.
4. A modified version of the reference-free aligner [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) is used to combine the alignments into base-level pangenome graphs that contain variants of all sizes
5. Chromosomal graphs are then combined and post-processed to reduce path complexity by collapsing redundant sequences. 

## 1.1 Fasta input files



## 1.2 bTaeGut.seqfile input file 

The main input file for the MC pipeline is the [SeqFile.txt]() file. Similarly to Progressive Cactus, this file contains the paths to all the fasta files we want to include in the pangenome. The main difference with the Progressive Cactus imput file is the absence of a guide three (Newick format), which is usually the first row in the Cactus file.

The file is organized like this:

| GenomeID.hap  | path/to/fasta |
| ------------- |:-------------:|
| bTaeGut1.1    | ~/2_fasta_files/..fasta |
| bTaeGut1.2    | ~/2_fasta_files/..fasta |
....etc

For the backbone reference, the pipeline currently doesn't support the presence of both haplotypes in the form bTaeGut1.1 and bTaeGut1.2. The coordinate space can only be based on a single haplotype and the variants inside the pangenome will be always referenced to the backbone reference, therefore it would not be possible to retrieve a VCF file with diploid genotype calls (e.g. 0/1) since one of the two haplotype it is indeed the reference. However, the pipeline will work using the IDs **bTaeGut1_hap1** and **bTaeGut1_hap2**. These will be considered as separate samples, but it's always convenient to include the alternate haplotype of the main reference as you will retain information about their variability for downstream analysis (e.g. reconstruction on an extinct species genome).

In addition to the chosen reference, one may specify additional assemblies with coordinates that can serve as a reference for graph decomposition. 

# References

1. Hickey, Glenn, et al. "Pangenome graph construction from genome alignments with Minigraph-Cactus." Nature biotechnology 42.4 (2024): 663-673.
