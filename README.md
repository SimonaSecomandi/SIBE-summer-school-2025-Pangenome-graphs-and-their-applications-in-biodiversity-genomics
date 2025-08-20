# Pangenome graphs and their application to biodiversity genomics

This repository contains the script and files for the pratical part of the course "Pangenome graphs and their application to biodiversity genomics" held the 8th of September 2025 in Ferrara, Italy, as part of the [SIBE summer school](https://sites.google.com/view/sibesummerschool/home-page).

The title of this course recalls our recent pangenomics review:
[Secomandi, S., Gallo, G., et al. Pangenome graphs and their applications in biodiversity genomics. Nat. Genet. 57, 13–26 (2025)](https://www.nature.com/articles/s41588-024-02029-6)

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

#If you plan to use the reference data, read the [README.md](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/reference_data/2_bTaeGut_pangenome/README.md) in the folder  [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/2_bTaeGut_pangenome) and follow the instructions.

### 0.2 Tools

Almost all the tools needed for this course are inside a conda environment.
Remember to **activate the conda environment** before running any command:

```
conda activate pangenomics
```
Other tools can be found here: 

## 1. Pangenome construction

To construct the pangenome we will use the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)<sup>2</sup>. The pipeline can be run as a whole or run step by step. 

#### Main steps of the MC pipeline<sup>1</sup>:

1. A user-selected reference genome is used as the initial backbone
2. The reference is progressively augmented with structural variation from the other genomes by minigraph<sup>3</sup>, a sequence-to-graph aligner, as a graph constructor. The resulting graph is SV only (>50 bp) (see **panel c** below)
3. All assemblies are aligned back to the graph with a minimap2-like algorithm that generates base-level alignments for each reference chromosome separately  (see **panel d** below)
4. The initian graph is split by chromosome and a modified version of the reference-free aligner Progressive Cactus<sup>4</sup> is used to combine the alignments into base-level pangenome graphs that contain variants of all sizes  (see **panel e** below)
6. Chromosomal graphs are then combined and post-processed to reduce path complexity by removing ("clipping") unaligned sequences such as the centromeres, leaving onlt the backbone reference in those regions and the other alleles are split accordingly (see **panel f** below)

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/figures/Figure1_Hickey_G_et_al.jpg" alt="drawing" width="600"/> <br/>
*Credits: [Hickey at al 2023](https://www.nature.com/articles/s41587-023-01793-w)*<br />

### 1.1 Fasta input files

We will create a pangenome for **chromosome 12** of two different **Zebra finch (*Taeniopygia guttata*)** individuals publicly available on NCBI.

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/figures/Zebra%20Finches_Michael_Lawton_Flickr.jpg" alt="drawing" width="300"/> <br/>
*Credits: [Flickr/Michael Lawton](https://www.flickr.com/photos/michaellawton/5712718319)*<br />

The **backbone reference** will be the new T2T reference genome [bTaeGut7.mat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048771995.1/), specifically chromosome 22 of maternal haplotype. We will also include chromosome 22 of paternal haplotype [bTaeGut7.pat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048772025.1/) and another individual's chromosome 22: [bTaeGut2.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051427915.1/) and [bTaeGut2.hap2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051428105.1/).

Differently from ProgressiveCactus, **MC ignores soft-masked bases**, so there is no need to repeat mask the genomes. For both, hard-masking is not recommended.

### 1.2 bTaeGut.seqfile input file 

The main input file for the MC pipeline is the [bTaeGut.seqfile ]() file. Similarly to Progressive Cactus, this text file contains user defined genome IDs and paths to the corresponding fasta file, one assembly per line. The main difference with the Progressive Cactus imput file is the absence of a guide three (Newick format), which is usually the first row in the Cactus file. 

The file is organized like this:

| GenomeID.hap  | path/to/fasta |
| ------------- |:-------------:|
| bTaeGut7_mat | 1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta |
| bTaeGut7_pat | 1_fasta_files/bTaeGut7_pat_chr22_CM109772.1.fasta |
| bTaeGut2.1 | 1_fasta_files/bTaeGut2_hap1_chr22_CM121079.1.fasta |
| bTaeGut2.2 | 1_fasta_files/bTaeGut2_hap2_chr22_CM121121.1.fasta |

**IMPORTANT:**
* Divergent haplotypes from the same individuals must be defined with ".1" and ".2" after the genome ID.
* The **backbone reference** (we will use *bTaeGut7.hap1*) can only be a single haplotype and **IT NEEDS TO BE CHROMOSOME-LEVEL**. Its chromosomes will be used as a template to align the other's genomes chromosomes/scaffolds and generate chromosomes graphs that will be then merged. This genome will be the main reference for the output VCF (containing the variants inside the pangenome) and will not appear as sample in the VCF. It will be also used as reference for downsteam analyses. 
* The pipeline currently does not support the presence of both haplotypes of the backbone individual in the form bTaeGut7.1 and bTaeGut7.2. The coordinate space for chromosome compartimentalization can only be based on a single haplotype (the backbone genome). However, the pipeline will work using the IDs **bTaeGut1_mat** (backbone ref) and **bTaeGut1_pat** (divergent haplotype of the backbone genome). You can also use "_hap1" and "_hap2", or what you prefer. These will be considered as separate samples, but it's always convenient to include the alternate haplotype of the main reference as you will retain information about their variability for downstream analyses (e.g. read mapping and variant calling). In the output VCF file, bTaeGut1_pat will appear as haploid (e.g. 0 instead of 1/0).

* In addition to the backbone reference, one may **specify additional assemblies as reference** when running the MC pipeline. The graph will be referenced to the backbone genome, but the additional genome's paths will serve as a reference for graph decomposition, for example, i.e. the pipeline will generate multiple VCF files, one referenced to the backbone reference and the others to the additional references (see below). The chromosome compartimentalization will still rely on the backbone references's chromosomes.

### 1.3 Running the MC pipeline

*The pipeline will take ~13 minutes to finish using 8 threads and 8 GB RAM.*

Activate the Cactus python environment:

```source /path/to/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate```

in my case: source /lustre/fs5/vgl/scratch/ssecomandi/BIN/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate

Run the following command in the main directory:

```
cactus-pangenome \
2_bTaeGut_pangenome/jobStore \
1_fasta_files/bTaeGut.seqfile \
--outDir 2_bTaeGut_pangenome \
--outName bTaeGut_pangenome \
--logFile 2_bTaeGut_pangenome/bTaeGut.log \
--reference bTaeGut7_mat bTaeGut7_pat \
--refContigs chr22 \
--vcf \
--vcfReference bTaeGut7_mat bTaeGut7_pat \
--filter 1 \
--gfa clip filter \
--gbz filter \
--giraffe filter \
--xg clip filter \
--chrom-vg clip filter \
--odgi clip full \
--viz clip full
```

* ```cactus-pangenome```: calls the MC pipeline 
* ```jobStore```: a directory used by the workflow management system to store intermediate files and job metadata
* ```bTaeGut.seqfile```: the file containing genomes IDs and fasta paths
* ```--mgMemory```: memory for minigraph construction. The memory can be a problem especially if you are running the jobs on a SLURM cluster, but we will set it anyways to be sure. Regarding the cores, the pipeline will use all those available.
* ```--consMemory```: memory for each cactus-consolidated job
* ```--indexMemory```: memory for indexing

Multiple flags can be set:

* ```--outDir```: folder where all the outputs will be generated
* ```--outName```: output prefixes
* ```--logFile```: path where to store the log file
* ```--reference```: the reference you want to use as a backbone followed by the other you want to use as reference for the output VCF
* ```--refContigs```: the reference chromosomes you want MC to base the chromosomes compartimentalization on. In our case it's just ```chr22``` of bTaeGut7_hap1
* ```--otherContigs``` (optional): tells MC to also consider unassembled scaffolds when generating the graph. They will be grouped in the same subgraph before merging the chromosomes. The name is used-defined. **Here is commented since we only have one chromosome.**
* ```--vcf```: tells MC to generate VCFs files referenced to the ```--vcfReference``` genomes 

In the following tags you can specify on which type of graph you want to operate on. Three type of graphs can be generated from the MC pipeline: 
1. **full graph**: all sequenced are included (default MC output)
2. **clip graph**: sequences not aligned the the SV-only graph produced by Minigraph and other sequences are removed to have a "cleaner" graph. 
3. **filter graph**: used for ```vg giraffe```, contains only nodes trasversed by X number of haplotypes (10% in the HPRC human pangenome)

* ```--filter 1```: removes nodes covered by less that 1 haplotype. We will use 1 since the 10% of our 4 haplotypes would be 0.4 
* ```--gfa clip filter```: output the graphs in the text-based Graphical Fragment Assembly (GFA) format, in addition to Cactus's native alignment format (HAL). This format is typically the most compatible for exchanging graphs between ```vg``` and other pangenome tools
* ```--gbz filter```: generate the .gbz graph file for these type of graphs. The filter grah in .gbz format is needed for ```vg giraffe```.
* ```--giraffe filter```: generate ```vg giraffe``` indexes for the filter graph (default)
* ```--chrom-vg clip filter```: generate the chromosome graphs in the Variation Graph (VG, .vg) format, usefull to use with the ```vg toolkit```. 
* ```--xg clip filter```: generate the .xg index file for the whole graph. It's a compressed, indexed representation of a variation graph, specifically optimized for fast path and graph traversal operations
*  ```--og full```: generate the graph in the ```odgi``` format  ```.og```. Usefull when running operations using ```odgi```
* ```--viz clip full```: generate an ```odgi viz``` 1D .png for each chromosome. ```odgi``` works better with ```full``` graphs, the presence off all sequences doesn't hinder the visualization. We will generate a .png also for the ```clip``` graph as a comparison. 

If you forget to set some of these flags don't worry, you can always generate ```.vg```, ```.og```, and other files later on.

### 1.4 The outputs

#### Graphs

```bTaeGut_pangenome.sv.gfa.gz```: SV-only graph output by Minigraph<br/>
```bTaeGut_pangenome.sv.gfa.fa.gz```: SVs included in the above Minigraph  graph in fasta format<br/>
```bTaeGut_pangenome.full.hal```: Cactus's native alignment format, can be used to convert to MAF (multiple Alignment Files) and use it for liftovers, create tracks on the UCSC browser or perfom conservation analyses, among others<br/>
```bTaeGut_pangenome.gfa.gz```: default ```clip``` graph in the default GFA format<br/>
```bTaeGut_pangenome.xg```: .xg index for the default ```clip``` graph<br/>
```bTaeGut_pangenome.og```: default ```clip``` graph in the ```odgi``` format<br/>
```bTaeGut_pangenome.full.og```: ```full``` graph in the ```odgi``` format, used for graph visualization with ```odgi```<br/>

```bTaeGut_pangenome.d1.gfa.gz```: ```filter``` graph in the default GFA format (1 = kept only nodes covered by at least 1 path)<br/>
```bTaeGut_pangenome.d1.gbz```: ```filter``` graph in GBZ format for ```vg giraffe```<br/>
```bTaeGut_pangenome.d1.xg```: .xg index for the ```filter``` graph<br/>
```bTaeGut_pangenome.d1.dist```: snarl distance index needed by ```vg giraffe```<br/>
```bTaeGut_pangenome.d1.min```: minimizer index needed by ```vg giraffe```<br/>
```bTaeGut_pangenome.d1.snarls```: end and start nodes for each bubble and nesting informations. Used by ```vg decontruct``` to generate the VCF files<br/>

#### VCFs

```bTaeGut_pangenome.raw.vcf.gz``` : main VCF file referenced to bTaeGut7_mat before normalization<br/>
```bTaeGut_pangenome.raw.vcf.gz.tbi```: index for the main VCF file referenced o bTaeGut7_mat before normalization<br/>
```bTaeGut_pangenome.vcf.gz```: final main VCF file referenced to bTaeGut7_mat)<br/>
```bTaeGut_pangenome.vcf.gz.tbi```: index for the final main VCF file referenced to bTaeGut7_mat<br/>
```bTaeGut_pangenome.bTaeGut7_pat.raw.vcf.gz```: VCF file referenced to bTaeGut7_pat before normalization<br/>
```bTaeGut_pangenome.bTaeGut7_pat.raw.vcf.gz.tbi```: index for VCF file referenced to bTaeGut7_pat before normalization<br/>
```bTaeGut_pangenome.bTaeGut7_pat.vcf.gz```: VCF file referenced to bTaeGut7_pat<br/>
```bTaeGut_pangenome.bTaeGut7_pat.vcf.gz.tbi```: index for VCF file referenced to bTaeGut7_pat<br/>

#### Additional files and folders

```bTaeGut.log```: log file of the MC run<br/>
```bTaeGut.seqfile```: copy of the input file<br/>
```bTaeGut_pangenome.stats.tgz```: clipping stats<br/>
```./bTaeGut_pangenome.chroms```: contains any graph type we speficied in the input for each chromosome<br/>
```./bTaeGut_pangenome.viz```: contains .png files generated with ```odgi viz```<br/>
```./chrom-alignments```: contains intermediate files<br/>
```./chrom-subproblems```: contains intermediate files<br/>

## 2. Pangenome evaluation

### 2.1 Statistics 
After the generation of the pangenome, the first thing to do is to check the statistics. This can be done with ```odgi stats``` or ```vg stats```.

Activate the conda environment with all our commands:

```conda activate SIBE_course```

We will generate general statistics for the pangenome using ```odgi```, starting from the ```.og``` ```clipped``` graph. As a reminder, this graph is the default MC graph and is a subgraph of the ```full``` graph in which sequences bigger than 10kb that were not aligned to the Minigraph SV-only graph and nodes that doesn't have edges on each side are removed. <br />

* **If you don't have an .og file** you can generate it as follows:

```bgzip -d -@ 8 2_bTaeGut_pangenome/bTaeGut_pangenome.gfa.gz```  <br />
```odgi build -t 8 -g 2_bTaeGut_pangenome/bTaeGut_pangenome.gfa -o 2_bTaeGut_pangenome/bTaeGut_pangenome.og```

* However, we asked MC to generate an odgi file with ```--og full clip```, so **we can directly generate statistics** for the ```clipped``` graph:

```odgi stats -S -i 2_bTaeGut_pangenome/bTaeGut_pangenome.og > 3_stats_and_viz/bTaeGut_pangenome.og.stats```

**The output:**

You can ```cat``` the output from ```odgi stats``` and look at the content:

```cat 3_stats_and_viz/bTaeGut_pangenome.og.stats```

| length|nodes|edges|paths|steps|
| ----- |:----:|:----:|:----:|:----:|
| 5904335|491480|670290|9|1198027

The graph has:
- **Length of 5.9 Mbp:** this is the pangenome graph sequences length, the sum of the lengths of all nodes in the graph &rarr; if you concatenate all node sequences, you would get ~5.9 Mbp of sequence
- **491k nodes**: the building blocks of the graph &rarr; each node corresponds to a contiguous DNA sequence
- **670k edges**:the connections between the nodes &rarr; edges describe how sequences can traverse from one node to another
- **9 paths**: a “path” usually corresponds to a genome/assembly/chromosome provided as input. See below for the explanation on why there are 9 paths instead of 4. 
- **1.19 M steps**: the total number of node visits across all paths

___

####  QUESTION 1: *is the pangenome graph bigger than the original reference sequence?* 

Generate the fasta index for the reference genome:

```samtools faidx 1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta```

Look at the .fai index file:

```cat 1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta.fai```

| chr | size   | offset  | linebases | linewidth |
| ----- |:----:|:----:|:----:|:----:|
chr22	| 5052704 |	7 | 70 | 71

In this index file, the chromosome size is in the second column.

* Original size of the backbone chromosome: **5.05 Mbp**
* Pangenome size: **5.9 Mbp**

#### ANSWER:  The pangenome size is bigger than the original reference. 
______

####  QUESTION 2: *of how much the reference was augmented by the other sequences?* 

*The size of a pangenome graph depends on the genome size of the respective species but is bound to be larger, as it incorporates accessory sequences from other individuals, and it is also influenced by the number and diversity of the individuals contributing to the pangenome as well as by the construction pipeline*<sup>1</sup> 

#### ANSWER: The other chromosomes augmented the reference by 0.9 Mbp
___

####  QUESTION 3: *Why do we have 4 input sequences but 9 paths? Shouldn’t there be only 4 paths?* 

Let's look at the paths inside the graph. You can list the paths with ```odgi paths```:

```odgi paths -L -i 2_bTaeGut_pangenome/bTaeGut_pangenome.og > 3_stats_and_viz/bTaeGut_pangenome.og.paths```

**Output:**

```cat 3_stats_and_viz/bTaeGut_pangenome.og.paths```

```
bTaeGut2#1#chr22#0[5771-877768]
bTaeGut2#1#chr22#0[889321-1697614]
bTaeGut2#1#chr22#0[1871689-2456870]
bTaeGut2#1#chr22#0[2530823-3686875]
bTaeGut2#1#chr22#0[3705235-3705399]
bTaeGut2#1#chr22#0[3743670-5050960]
bTaeGut2#2#chr22#0[0-4660045]
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22[6348-4838699]
```

The format is: ```sample#hap#chrom#coords[start-end]```

* *bTaeGut7_mat#0#chr22* = **Path 1** &rarr; the backbone reference is in one piece (used to map the other chromosomes) 
* bTaeGut7_pat#0#chr22[6348-4838699] = **Path 2** &rarr; the alternate haplotype of the backbone reference is aligning in a big piece with start coordinate at 6348 bp and end coordinate at 4838699 bp
* bTaeGut2#1#chr22#0[*] = **Paths 3-8** &rarr; haplotype 1 of the second individual is not aligning contigously and it's split on 6 different pieces
* bTaeGut2#2#chr22#0[0-4660045] = **Path 9** &rarr; haplotype 2 of the second individual is also aligning in a big piece

#### ANSWER: if the sequences doesn't align contiguously to the backbone reference (e.g. in the presence of repetitive regions) the other chromosomes might be split in multiple segments and counted as different paths, 9 in our case. Functionally, all bTaeGut2#1#chr22#0[*] lines together represent one biological chromosome, but they are split across multiple path segments.
___

#### QUESTION 4: *Does this occur in the full graph too, or only after clipping?*

Generate the stats for the ```full``` graph (we already have a ```.full.og``` file, we asked MC to generate it with ```--og full clip```):

```odgi stats -S -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og > 3_stats_and_viz/bTaeGut_pangenome.full.og.stats```

```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.stats```

| length|nodes|edges|paths|steps|
| ----- |:----:|:----:|:----:|:----:|
| 6936749|492938|672043|4|1199841

Let's look at the paths inside the graph. You can list the paths with ```odgi paths```:

```odgi paths -L -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og > 3_stats_and_viz/bTaeGut_pangenome.full.og.paths```

**Output:**

```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.paths```

```
bTaeGut2#1#chr22#0
bTaeGut2#2#chr22#0
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22
```

#### ANSWER: in the full graph we have exactly 4 paths and the sequence length is way bigger than the clipped (6.9 Mbp). 

* **Full graph**: keeps each input chromosome/haplptype as a single continuous path, 4 in total
* **Clip graph**: removed sequences that don't overlap across genomes and long chromosomes get split into multiple segments and that's why bTaeGut2#1#chr22#0 appears in 7 separate path segments. The overall graph length is also smaller than the full (~6.9 Mbp &rarr; ~5.9 Mbp))

___

#### QUESTION 5: *why is clipping necessary?*

Clipping is necessary to reduce the complexity of the graph and make downstream analyses more feasible. For example, the very repetitive centromeric regions are clipped away and this is necessary to avoid the generation of complex loops that will hinder downstream analyses.

### 2.2 Subsampling and visualization 

Visualization is important to get an idea of the structure of the graph and the inspection of homology relationships and variation between the genomes, providing insights on the latent biological data.

#### odgi viz<sup>5</sup>

In the MC command we specified to generate ```odgi viz``` graphs and these can be found in the folder XX. A .png file was generated for each chromosome (one in our particular case) starting from the ```full graph``` and the ```clip graph```. 

Let's look at the ```full graph```'s png. Each line represents a different chromosome' path with their genome.ID on the right side. Each path is coloured when passing through a node and edges representing variation are represented by black lines in the bottom of the figure.

There is the possibility of ordering the paths as we want, e.g. to have the backbone refreence on the top.

To do this we can use the ```3_stats_and_viz/bTaeGut_pangenome.full.og.paths``` file, reordering the paths as we prefer. 


```
(grep "bTaeGut7_mat#0#chr22" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut7_pat#0#chr22" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut2#1#chr22#0" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut2#2#chr22#0" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths) \
> 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths
```
Let's look at the new path file to check if it's the correct order:

```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths```

```
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22
bTaeGut2#1#chr22#0
bTaeGut2#2#chr22#0
```

Now we can generate a new ```odgi viz``` figure with the ordered paths:

```odgi viz -x 1500 -y 500 -a 10 -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og -p 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths -o 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.viz.png```

#### SequenceTubeMap<sup>6</sup>

Another useful tool for visualizing pangenome graphs is SequenceTubeMap<sup>6</sup>. It visualizes the graph in ```.vg``` format using the same linear visualization as ```odgi viz```, but variability among genomes is displayed differently and it can be inspected interactively.

Unfortunately, the online demo doesn't work right now. To visualize a specific vg file without uploading it on the webpage, it is possible to launch a server which provides the data to SequenceTubeMap<sup>6</sup>. See instructions on the [SequenceTubeMap GitHub page](https://github.com/vgteam/sequenceTubeMap). In this course we will just prepare the files and you can try at home. 

First, we will chunk the graph in a smaller piece to be able to visualize it fast.

1. Chunk the graph: 

```vg chunk -t 8 -c 1 -x 2_bTaeGut_pangenome/bTaeGut_pangenome.xg -p bTaeGut7_mat#0#chr22:0-100000 -O vg > 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.vg```

*Ignore the warning: "warning[vg chunk]: the vg-protobuf format is DEPRECATED. you probably want to use PackedGraph (pg) instead"*

* ```-c, --context-steps N```: expand the context of the chunk this many node steps [1]

2. Index the new .vg chunk:

```vg convert -t 8 -x 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.vg > 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.xg```

2. Upload the files in the [online demo](https://vgteam.github.io/sequenceTubeMap/). You can find the files in this github repository, download it directly from her the the computer. Follow the intructions:
* Go to the [sequenceTubeMap demo page](https://vgteam.github.io/sequenceTubeMap/). 
* Select "Custom" from the "Data" drop down menu > Click on "Configure Tracks" > click the "+" button > leave "graph" but change the "mounted" with "upload" > select the file from the Download folder > close 

Inspect the graph:
- How many variants do you see?
- ....

Unfortunately, the online demo does not work with custom data. Let's look at the examples already present in the online demo. While I'll show you our pangenome by sharing my screen.


## 3. Pangenome-embedded variants

We looked at the structure of the graph and the variants inside the pangenome, but *how can I look at them in a canonical way and use them for downstream analysis?*

The MC pipeline produces **VCF files** referenced to the backbone reference and other genomes you specified in the command. You can find information about the VCF format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

The process through which the variants are defined is called **graph decomposition**, the process of breaking down a pangenome graph into smaller, more manageable subgraphs or components (snarls or bubbles). 

You might have noticed "raw" VCF files. These are those directly outputted by vg deconstruct inside that are then normalized and postprocessed automatically.

Let's look at the VCF referenced to our backbone reference.

___

#### QUESTION 6: *how many samples do you see in the VCF?*

You can find information about the samples in the last row of the header file (i.e. the last line starting with "#") after the column "FORMAT".

Run the following command in the main directory:

```bcftools view 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | head -20```

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CONFLICT,Number=.,Type=String,Description="Sample names for which there are multiple paths in the graph with conflicting alleles">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=chr22,length=5052704>
##bcftools_viewVersion=1.21+htslib-1.21
##bcftools_viewCommand=view 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz; Date=Wed Aug 20 15:34:06 2025
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	bTaeGut2	bTaeGut7_pat
chr22	14063	>374>377	C	G	60	.	AC=1;AF=0.5;AN=2;AT=>374>375>377,>374>376>377;NS=2;LV=0	GT	1|.	0
chr22	14268	>377>379	T	TT	60	.	AC=1;AF=0.5;AN=2;AT=>377>379,>377>378>379;NS=2;LV=0	GT	1|.	0
chr22	14304	>379>382	G	A	60	.	AC=1;AF=0.5;AN=2;AT=>379>380>382,>379>381>382;NS=2;LV=0	GT	1|.	0
chr22	14312	>382>385	C	A	60	.	AC=1;AF=0.5;AN=2;AT=>382>383>385,>382>384>385;NS=2;LV=0	GT	1|.	0
chr22	14316	>385>388	TA	CT	60	.	AC=1;AF=0.5;AN=2;AT=>385>386>388,>385>387>388;NS=2;LV=0	GT	1|.	0
```

#### ANSWER: there are just two samples, one is the haploid individual bTaeGut2, with diploid genotype calls (i.e. 1|1), the other is the paternal haplotype of the backbone reference, which is indeed haploid (i.e. 0)

___

#### QUESTION 7: *why the diploid individual's genotypes are separed by a "|"?*

Usually, when you call the variants from short-read mapping, for example, you have calls like this one "1\1".

in a VCF:
* ```/``` &rarr; **unphased** = we know the two alleles are present but don’t know which belongs to which haplotype
* ```|``` &rarr; **phased** = we know which allele is on haplotype 1 vs haplotype 2

#### ANSWER: the MC pipeline can output phased genotype because when you build a pangenome with chromosome-level haplotypes, each haplotype is represented explicitly as a path through the graph. The variant calling is performed by decomposing the graph into “bubbles” (alternative allelic paths) and since the graph keeps track of which haplotype path traverses which node in the bubble, the software can directly assign alleles to haplotypes.
____

# References

1. Secomandi, Simona, et al. "Pangenome graphs and their applications in biodiversity genomics." Nature Genetics 57.1 (2025): 13-26.
2. Hickey, Glenn, et al. "Pangenome graph construction from genome alignments with Minigraph-Cactus." Nature biotechnology 42.4 (2024): 663-673.
3. Li, Heng, Xiaowen Feng, and Chong Chu. "The design and construction of reference pangenome graphs with minigraph." Genome biology 21.1 (2020): 265.
4. Armstrong, Joel, et al. "Progressive Cactus is a multiple-genome aligner for the thousand-genome era." Nature 587.7833 (2020): 246-251.
5. Guarracino, Andrea, et al. "ODGI: understanding pangenome graphs." Bioinformatics 38.13 (2022): 3319-3326.
6. Beyer, Wolfgang, et al. "Sequence tube maps: making graph genomes intuitive to commuters." Bioinformatics 35.24 (2019): 5318-5320.
