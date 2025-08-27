# Pangenome graphs and their application to biodiversity genomics

This repository contains the script and files for the pratical part of the course "Pangenome graphs and their application to biodiversity genomics" held the 8th of September 2025 in Ferrara, Italy, as part of the [SIBE summer school](https://sites.google.com/view/sibesummerschool/home-page).

The title of this course recalls our recent pangenomics review:
[Secomandi, S., Gallo, G., et al. Pangenome graphs and their applications in biodiversity genomics. Nat. Genet. 57, 13–26 (2025)](https://www.nature.com/articles/s41588-024-02029-6)

## Introduction

This course will...

## Table of content

**0. First steps**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.1 Files and folders<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.2 Tools<br />
**1. Pangenome construction**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1 Fasta input files<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2 bTaeGut.seqfile input file<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3 Running the MC pipeline<br />
**2. Pangenome evaluation and statistics**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.1 Generate graph statistics<br />
**3. Visualization and subsampling**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.1 odgi viz<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2 SequenceTubeMap<br />
**4. Pangenome-embedded small variants**<br />
**5. Mapping of short-reads data with vg giraffe**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.1 Align the reads with vg giraffe<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.2 Project graph alignments onto a linear reference with vg surject<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5.3 Visualize aligned reads with SequenceTubeMap6<br />

## 0. First steps

### 0.1 Files and folders

The main directory of this repository contains all the input files needed for the exercises, as well as empty folders where your outputs will be written. All commands can be run from the main directory and they will output the data in the correct folders automatically.

A separate folder, [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data), contains all the expected outputs organized using the same folder structure. You can use these reference files to compare results with the expected outcome, complete exercises after the session at home and troubleshoot and debug each step.

#### 0.1.1 Copy this repository in your folder:

```
git clone https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics.git
```

#### 0.1.2 Download "big" files

Unfortunately GitHub doesn't allow the upload of files > 100 Gbp. You can download these files from GoogleDrive at this [link](https://drive.google.com/drive/folders/1UQXGIKwR2ErMcv6pT097QnPD9XzjChrL?usp=sharing).

Download them and place them in the correct folder ```reference_data/5.1_vg_giraffe``` and ```reference_data/5.2_bwa_mem```.

### 0.2 Tools

Almost all the tools needed for this course are inside a conda environment, except ```Minigraph-Cactus```.

Remember to **activate the Cactus environment** before running ```Minigraph-Cactus```:

```source /path/to/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate```

Remember to **activate the conda environment** before running any command after the pangenome contrustion:

```conda activate SIBE_course```

You will find these commands throught the excercises when needed. 

**Rerember, all the commands MUST be run from the main directory!**

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

We will create a pangenome for a chunk of **chromosome 22** (the first Mbp) of two different **Zebra finch (*Taeniopygia guttata*)** individuals publicly available on NCBI.

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/figures/Zebra%20Finches_Michael_Lawton_Flickr.jpg" alt="drawing" width="300"/> <br/>
*Credits: [Flickr/Michael Lawton](https://www.flickr.com/photos/michaellawton/5712718319)*<br />

The **backbone reference** will be the new T2T reference genome [bTaeGut7.mat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048771995.1/), specifically chromosome 22 of maternal haplotype. We will also include chromosome 22 of paternal haplotype [bTaeGut7.pat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048772025.1/) and another individual's chromosome 22: [bTaeGut2.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051427915.1/) and [bTaeGut2.hap2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051428105.1/).

Differently from ProgressiveCactus, **MC ignores soft-masked bases**, so there is no need to repeat mask the genomes. For both, hard-masking is not recommended.

### 1.2 bTaeGut.seqfile input file 

The main input file for the MC pipeline is the [bTaeGut.seqfile](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/1_fasta_files/bTaeGut.seqfile) file in the folder ```1_fasta_files/```. Similarly to ```Progressive Cactus```, this text file contains user defined genome IDs and paths to the corresponding fasta file, one assembly per line. The main difference with the Progressive Cactus imput file is the absence of a guide three (```Newick``` format), which is usually the first row in the ```Cactus``` input file. 

The file is organized like this:

| GenomeID.hap  | path/to/fasta |
| ------------- |:-------------:|
| bTaeGut7_mat | 1_fasta_files/bTaeGut7_mat_chr22_1Mbp_NC_133047.1.fasta |
| bTaeGut7_pat | 1_fasta_files/bTaeGut7_pat_chr22_1Mbp_CM109772.1.fasta |
| bTaeGut2.1 | 1_fasta_files/bTaeGut2_hap1_chr22_1Mbp_CM121079.1.fasta |
| bTaeGut2.2 | 1_fasta_files/bTaeGut2_hap2_chr22_1Mbp_CM121121.1.fasta |

**IMPORTANT:**
* Divergent haplotypes from the same individuals must be defined with ".1" and ".2" after the genome ID.
* The **backbone reference** (we will use *bTaeGut7.hap1*) can only be a single haplotype and **IT NEEDS TO BE CHROMOSOME-LEVEL**. Its chromosomes will be used as a template to align the other's genomes chromosomes/scaffolds and generate chromosomes graphs that will be then merged. This genome will be the main reference for the output VCF (containing the variants inside the pangenome) and will not appear as sample in the VCF. It will be also used as reference for downsteam analyses. 
* The pipeline currently does not support the presence of both haplotypes of the backbone individual in the form bTaeGut7.1 and bTaeGut7.2. The coordinate space for chromosome compartimentalization can only be based on a single haplotype (the backbone genome). However, the pipeline will work using the IDs **bTaeGut1_mat** (backbone ref) and **bTaeGut1_pat** (divergent haplotype of the backbone genome). You can also use "_hap1" and "_hap2", or what you prefer. These will be considered as separate samples, but it's always convenient to include the alternate haplotype of the main reference as you will retain information about their variability for downstream analyses (e.g. read mapping and variant calling). In the output VCF file, bTaeGut1_pat will appear as haploid (e.g. 0 instead of 1/0).

* In addition to the backbone reference, one may **specify additional assemblies as reference** when running the MC pipeline. The graph will be referenced to the backbone genome, but the additional genome's paths will serve as a reference for graph decomposition, for example, i.e. the pipeline will generate multiple VCF files, one referenced to the backbone reference and the others to the additional references (see below). The chromosome compartimentalization will still rely on the backbone references's chromosomes.

### 1.3 Running the MC pipeline

*The pipeline will take ~13 minutes to finish using 8 threads and 8 GB RAM.*

#### 1.3.1 Activate the Cactus python environment:

**RUN:**
```source /path/to/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate```

in my case: source /lustre/fs5/vgl/scratch/ssecomandi/BIN/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate

#### 1.3.2 Run the pipeline:
   
**RUN**:

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

The pipeline will take around 3 minutes to run. You can find the log (the same you saw in stdout) here:  
Copy and paste the outputs from the [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/2_bTaeGut_pangenome) folder if there is no time to run the command:
```
cp -r reference_data/2_bTaeGut_pangenome/* 2_bTaeGut_pangenome
```

**The inputs:**

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

**The outputs:**

Graphs:

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

VCFs:

```bTaeGut_pangenome.raw.vcf.gz``` : main VCF file referenced to bTaeGut7_mat before normalization<br/>
```bTaeGut_pangenome.raw.vcf.gz.tbi```: index for the main VCF file referenced o bTaeGut7_mat before normalization<br/>
```bTaeGut_pangenome.vcf.gz```: final main VCF file referenced to bTaeGut7_mat)<br/>
```bTaeGut_pangenome.vcf.gz.tbi```: index for the final main VCF file referenced to bTaeGut7_mat<br/>
```bTaeGut_pangenome.bTaeGut7_pat.raw.vcf.gz```: VCF file referenced to bTaeGut7_pat before normalization<br/>
```bTaeGut_pangenome.bTaeGut7_pat.raw.vcf.gz.tbi```: index for VCF file referenced to bTaeGut7_pat before normalization<br/>
```bTaeGut_pangenome.bTaeGut7_pat.vcf.gz```: VCF file referenced to bTaeGut7_pat<br/>
```bTaeGut_pangenome.bTaeGut7_pat.vcf.gz.tbi```: index for VCF file referenced to bTaeGut7_pat<br/>

Additional files and folders:

```bTaeGut.log```: log file of the MC run<br/>
```bTaeGut.seqfile```: copy of the input file<br/>
```bTaeGut_pangenome.stats.tgz```: clipping stats<br/>
```./bTaeGut_pangenome.chroms```: contains any graph type we speficied in the input for each chromosome<br/>
```./bTaeGut_pangenome.viz```: contains .png files generated with ```odgi viz```<br/>
```./chrom-alignments```: contains intermediate files<br/>
```./chrom-subproblems```: contains intermediate files<br/>

## 2. Pangenome evaluation and statistics

After the generation of the pangenome, the first thing to do is to check the statistics. This can be done with ```odgi stats``` or ```vg stats```.

##### Deactivate the Cactus environment:

**RUN:**
```deactivate```

##### Activate the conda environment with all our commands:

**RUN:**
```conda activate SIBE_course```

#### 2.1 Generate graph statistics

We will generate general statistics for the pangenome using ```odgi```, starting from the ```.og``` ```clipped``` graph. As a reminder, this graph is the default MC graph and is a subgraph of the ```full``` graph in which sequences bigger than 10kb that were not aligned to the Minigraph SV-only graph and nodes that doesn't have edges on each side are removed. <br />

**If you don't have an .og file** you can generate it as follows

```bgzip -d -@ 8 2_bTaeGut_pangenome/bTaeGut_pangenome.gfa.gz```  <br />
```odgi build -t 8 -g 2_bTaeGut_pangenome/bTaeGut_pangenome.gfa -o 2_bTaeGut_pangenome/bTaeGut_pangenome.og```

However, we asked MC to generate an odgi file with ```--og full clip```, so **we can directly generate statistics** for the ```clipped``` graph

**RUN:**
```odgi stats -S -i 2_bTaeGut_pangenome/bTaeGut_pangenome.og > 3_stats_and_viz/bTaeGut_pangenome.og.stats```

**The inputs:**

* ```-S```: print a summary of the graph properties and dimensions
* ```-i```: the input ```.og``` variantion graph

**The output:**

You can ```cat``` the output from ```odgi stats``` and look at the content.

**RUN:**
```cat 3_stats_and_viz/bTaeGut_pangenome.og.stats```

| length|nodes|edges|paths|steps|
| ----- |:----:|:----:|:----:|:----:|
| 1123880|54958|74832|5|136900

The graph has:
- **Length ~1.12 Mbp:** this is the pangenome graph sequences length, the sum of the lengths of all nodes in the graph &rarr; if you concatenate all node sequences, you would get ~1.12 Mbp of sequence
- **~54.9k nodes**: the building blocks of the graph &rarr; each node corresponds to a contiguous DNA sequence
- **~670k edges**:the connections between the nodes &rarr; edges describe how sequences can traverse from one node to another
- **~5 paths**: a “path” usually corresponds to a genome/assembly/chromosome provided as input. See below for the explanation on why there are 5 paths instead of 4. 
- **~136k steps**: the total number of node visits across all paths

___

####  QUESTION 1: *is the pangenome graph bigger than the original reference sequence?* 

#### 1. Generate the fasta index for the reference genome:

**RUN:**
```samtools faidx 1_fasta_files/bTaeGut7_mat_chr22_1Mb_NC_133047.1.fasta```

#### 2. Look at the .fai index file:

**RUN:**
```cat 1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta.fai```

| chr | size   | offset  | linebases | linewidth |
| ----- |:----:|:----:|:----:|:----:|
chr22	| 1000000 |	7 | 1000000 | 1000001

In this index file, the chromosome size is in the second column.

* Original size of the backbone chromosome: **1.12 Mbp**
* Pangenome size: **1 Mbp**

#### ANSWER:  The pangenome size is bigger than the original reference. 

______

####  QUESTION 2: *by how much was the reference augmented by the other sequences??* 

*The size of a pangenome graph depends on the genome size of the respective species but is bound to be larger, as it incorporates accessory sequences from other individuals, and it is also influenced by the number and diversity of the individuals contributing to the pangenome as well as by the construction pipeline*<sup>1</sup> 

#### 1. calculate the difference between the pangenome length and the reference length 

**pangenome length - reference length** <br/>
1123880 bp - 1000000 bp = 123880 bp

#### ANSWER: The other chromosomes augmented the reference by ~123 Kbp

___

####  QUESTION 3: *Why do we have 4 input sequences but 5 paths?* 

#### 1. Let's look at the paths inside the graph. You can list the paths with ```odgi paths```.

**RUN:**
```odgi paths -L -i 2_bTaeGut_pangenome/bTaeGut_pangenome.og > 3_stats_and_viz/bTaeGut_pangenome.og.paths```

#### 2. Look at the path list:

**RUN:**
```cat 3_stats_and_viz/bTaeGut_pangenome.og.paths```

```
bTaeGut2#1#chr22#0[5771-877768]
bTaeGut2#1#chr22#0[889321-995863]
bTaeGut2#2#chr22#0
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22[6348-993008]
```

The naming format is: ```sample#hap#chrom#coords[start-end]```, also called [PanSN](https://github.com/pangenome/PanSN-spec) format. 

* *bTaeGut7_mat#0#chr22* = **Path 1** &rarr; the backbone reference is in one piece (used to map the other chromosomes) 
* bTaeGut7_pat#0#chr22[6348-993008] = **Path 2** &rarr; the alternate haplotype of the backbone reference is aligning in a big piece with start coordinate at 6348 bp and end coordinate at 993008 bp
* bTaeGut2#1#chr22#0[*] = **Paths 2-3** &rarr; haplotype 1 of the second individual is not aligning contigously and it's split into 2 different pieces
* bTaeGut2#2#chr22#0 = **Path 4** &rarr; haplotype 2 of the second individual is also aligning in a big piece

#### ANSWER: if the sequences doesn't align contiguously to the backbone reference (e.g. in the presence of repetitive regions) the other chromosomes might be split in multiple segments and counted as different paths, 5 in our case. Functionally, all bTaeGut2#1#chr22#0[*] lines together represent one biological chromosome, but they are split across multiple path segments.
___

#### QUESTION 4: *Does this occur in the full graph too, or is this a result of clipping?*

#### 1. Generate the stats for the ```full``` graph (we already have a ```.full.og``` file, we asked MC to generate it with ```--og full clip```).

**RUN:**
```odgi stats -S -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og > 3_stats_and_viz/bTaeGut_pangenome.full.og.stats```

**RUN:**
```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.stats```

| length|nodes|edges|paths|steps|
| ----- |:----:|:----:|:----:|:----:|
| 1150037|55483|75569|4|137676

#### 2. Let's now look at the paths inside the graph. You can list the paths with ```odgi paths```.

**RUN:**
```odgi paths -L -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og > 3_stats_and_viz/bTaeGut_pangenome.full.og.paths```

**RUN:**
```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.paths```

```
bTaeGut2#1#chr22#0
bTaeGut2#2#chr22#0
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22
```

#### ANSWER: in the full graph we have exactly 4 paths and the sequence length is way bigger than the clipped (6.9 Mbp). 

* **Full graph**: keeps each input chromosome/haplptype as a single continuous path, 4 in total
* **Clip graph**: removed sequences that don't overlap across genomes and long chromosomes get split into multiple segments and that's why bTaeGut2#1#chr22#0 appears in 2 separate path segments. The overall graph length is also smaller than the full (~1.12 Mbp &rarr; ~1.15 Mbp))

___

#### QUESTION 5: *why is clipping necessary?*

Clipping is necessary to reduce the complexity of the graph and make downstream analyses more feasible. For example, the very repetitive centromeric regions are clipped away and this is necessary to avoid the generation of complex loops that will hinder downstream analyses.

___

## 3. Visualization and subsampling

Visualization is important to get an idea of the structure of the graph and the inspection of homology relationships and variation between the genomes, providing insights on the latent biological data.

### 3.1 ```odgi viz```<sup>5</sup>

In the MC command we specified to generate ```odgi viz``` graphs and these can be found in the folder 2_bTaeGut_pangenome/bTaeGut_pangenome.viz. A .png file was generated for each chromosome (one in our particular case) starting from the ```full graph``` and the ```clip graph```. 

#### 3.1.1 Let's look at the ```full graph```'s png (```2_bTaeGut_pangenome/bTaeGut_pangenome.viz/chr22.full.viz.png```) 

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/reference_data/2_bTaeGut_pangenome/bTaeGut_pangenome.viz/chr22.full.viz.png" alt="drawing" width="1000"/> <br/>

Each line represents a different chromosome' path with their genome.ID on the right side. Each path is coloured when passing through a node and edges representing variation are represented by black lines in the bottom of the figure.

#### 3.1.2 Change the path order in the .png

There is the possibility of ordering the paths as we want, e.g. to have the backbone reference at the top.

To do this we can use the ```3_stats_and_viz/bTaeGut_pangenome.full.og.paths``` file, reordering the paths as we prefer. 

**RUN:**
```
(grep "bTaeGut7_mat#0#chr22" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut7_pat#0#chr22" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut2#1#chr22#0" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths; \
grep "bTaeGut2#2#chr22#0" 3_stats_and_viz/bTaeGut_pangenome.full.og.paths) \
> 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths
```

Let's look at the new path file to check if it's the correct order:

**RUN:**
```cat 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths```

```
bTaeGut7_mat#0#chr22
bTaeGut7_pat#0#chr22
bTaeGut2#1#chr22#0
bTaeGut2#2#chr22#0
```

The order is correct! Now we can generate a new ```odgi viz``` figure with the ordered paths:

**RUN:**
```odgi viz -x 1500 -y 500 -a 10 -i 2_bTaeGut_pangenome/bTaeGut_pangenome.full.og -p 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.paths -o 3_stats_and_viz/bTaeGut_pangenome.full.og.sort.viz.png```

The flags are to indicate the image width (```-x```) and height (```-y```).

Look at the new .png (```3_stats_and_viz/bTaeGut_pangenome.full.og.sort.viz.png```):

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/reference_data/3_stats_and_viz/bTaeGut_pangenome.full.og.sort.viz.png" alt="drawing" width="1000"/> <br/>

Now the paths are in the correct order.

### 3.2 ```SequenceTubeMap```<sup>6</sup>

Another useful tool for visualizing pangenome graphs is ```SequenceTubeMap```<sup>6</sup>. It visualizes the graph in ```.vg``` format using the same linear visualization as ```odgi viz```, but variability among genomes is displayed differently and it can be inspected interactively.

There is an [online demo](https://vgteam.github.io/sequenceTubeMap/) that can be used to upload small files. However, sometime it doesn't work (since it's just a demo) and most probably does not support multiple users loading data at the same time. Let's try it. If it doesn't work, you can look at my screen and the screenshots and videos I uploaded in the folder ```3_stats_and_viz/sequenceTubeMap/```:

To visualize a specific graph file without uploading it on the webpage, it is possible to launch a server which provides the data to ```SequenceTubeMap```<sup>6</sup>. See instructions on the [```SequenceTubeMap``` GitHub page](https://github.com/vgteam/```SequenceTubeMap```). In this course we will just prepare the files and you can try later at home.

First, we will chunk the graph in a smaller piece to be able to visualize it fast.

#### 3.2.1 Chunk the graph 

You can visualize a ```.vg``` graph and it's index ```.xg``` with ```SequenceTubeMap```. To aid visualization and avoid using too much memory, we will subsample the graph at specific coordinates.

**RUN:**
```vg chunk -t 4 -c 1 -x 2_bTaeGut_pangenome/bTaeGut_pangenome.xg -p bTaeGut7_mat#0#chr22:0-100000 -O vg > 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.vg```

The flag ```-c, --context-steps N``` tells ```vg chunk``` to expand the context of the chunk this many node steps.

*Ignore the warning: "warning[vg chunk]: the vg-protobuf format is DEPRECATED. you probably want to use PackedGraph (pg) instead"*

#### 3.2.2 Index the new ```.vg``` chunk

**RUN:**
```vg convert -t 8 -x 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.vg > 3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.xg```

#### 3.2.3 Look at the paths inside the chunk

We can use ```vg path``` since we are looking at the ```.xg``` file.

**RUN:**
```vg paths -L -x bTaeGut_pangenome.chunk.100Kb.xg```
```
bTaeGut7_pat#0#chr22[6348-38182]
bTaeGut7_pat#0#chr22[38183-91119]
bTaeGut7_mat#0#chr22[0-100017]
bTaeGut2#1#chr22#5771[13553-35791]
bTaeGut2#1#chr22#5771[42325-93155]
bTaeGut2#2#chr22#0[0-48874]
bTaeGut2#1#chr22#5771[0-2494]
bTaeGut2#1#chr22#5771[35792-39253]
bTaeGut2#1#chr22#5771[4542-13552]
```

As you can see we have more paths, this happens because the chunk may include disconnected pieces of a path (because nodes outside the requested range were dropped), so each connected run becomes its own [start-end] fragment. ```vg chunk``` splits each original path into the pieces that lie inside the subgraph we extracted. 

#### 3.2.4 Upload the files in the [online demo](https://vgteam.github.io/SequenceTubeMap/). 

You can find the files in this github repository, download it directly from her the the computer:
* ```3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.vg```
* ```3_stats_and_viz/bTaeGut_pangenome.chunk.100Kb.xg```

Follow the instructions:
* Go to the [SequenceTubeMap demo page](https://vgteam.github.io/SequenceTubeMap/). 
* Select "Custom" from the "Data" drop down menu > Click on "Configure Tracks" > click the "+" button > leave "graph" but change "mounted" with "upload" > select the .xg file from the Download folder > close using the "x" in the upper right corner

## 4. Pangenome-embedded small variants

We looked at the structure of the graph and the variants inside the pangenome, but *how can I look at them in a canonical way and use them for downstream analysis?*

The MC pipeline produces **VCF files** referenced to the backbone reference and other genomes you specified in the command. You can find information about the ```VCF ```format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

The process through which the variants are defined is called **graph decomposition**, the process of breaking down a pangenome graph into smaller, more manageable subgraphs or components (snarls or bubbles). 

You might have noticed "raw" ```VCF``` files. These are those directly outputted by ```vg deconstruct``` inside that are then normalized and postprocessed automatically.

Let's look at the ```VCF``` referenced to our backbone reference.

___

#### QUESTION 6: *how many samples do you see in the VCF?*

You can find information about the samples in the last row of the ```VCF``` ```header``` (i.e. the last line starting with "#") after the column ```"FORMAT"```.

**RUN:**
```bcftools view 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | head -20```

```
#fileformat=VCFv4.2
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
##contig=<ID=chr22,length=1000000>
##bcftools_viewVersion=1.21+htslib-1.21
##bcftools_viewCommand=view 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz; Date=Wed Aug 27 10:29:53 2025
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	bTaeGut2	bTaeGut7_pat
chr22	14063	>351>354	C	G	60	.	AC=1;AF=0.5;AN=2;AT=>351>352>354,>351>353>354;NS=2;LV=0	GT	1|.	0
chr22	14268	>354>356	T	TT	60	.	AC=1;AF=0.5;AN=2;AT=>354>356,>354>355>356;NS=2;LV=0	GT	1|.	0
chr22	14304	>356>359	G	A	60	.	AC=1;AF=0.5;AN=2;AT=>356>357>359,>356>358>359;NS=2;LV=0	GT	1|.	0
chr22	14312	>359>362	C	A	60	.	AC=1;AF=0.5;AN=2;AT=>359>360>362,>359>361>362;NS=2;LV=0	GT	1|.	0
chr22	14316	>362>365	TA	CT	60	.	AC=1;AF=0.5;AN=2;AT=>362>363>365,>362>364>365;NS=2;LV=0	GT	1|.	0
```

#### ANSWER: there are just two samples, one is the haploid individual bTaeGut2, with diploid genotype calls (i.e. ```1|1```), the other is the paternal haplotype of the backbone reference, which is indeed haploid (i.e. ```0```)

___

#### QUESTION 7: *why the diploid individual's genotypes are separed by a "|"?*

Usually, when you call the variants from short-read mapping, for example, you have calls like this one ```1\1```.

in a ```VCF```:
* ```/``` &rarr; **unphased** = we know the two alleles are present but don’t know which belongs to which haplotype
* ```|``` &rarr; **phased** = we know which allele is on haplotype 1 vs haplotype 2

#### ANSWER: the ```MC pipeline``` can output phased genotype because when you build a pangenome with chromosome-level haplotypes, each haplotype is represented explicitly as a path through the graph. The variant calling is performed by decomposing the graph into “bubbles” (alternative allelic paths) and since the graph keeps track of which haplotype path traverses which node in the bubble, the software can directly assign alleles to haplotypes.

____

#### QUESTION 8: *how many SNPs and how many INDELS our pangenome includes?*

Let's look at biallelic SNPs and INDELS, which are those suitable to generate Principal Component Analysis (PCA) plots and other population genomics analysis.

* To count biallelic SNPs **RUN:**

**RUN:**
```bcftools view -v snps --max-alleles 2 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | grep -v "^#" | wc -l```

* To count biallelic INDELs **RUN:**

```bcftools view -v indels --max-alleles 2 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | grep -v "^#" | wc -l```

* To count the number of insertions and deletions **RUN:**

```bcftools view -v indels --max-alleles 2 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | awk 'length($4) > length($5)' | wc -l```

* To count the number of deletions **RUN:**

```bcftools view -v indels --max-alleles 2 2_bTaeGut_pangenome/bTaeGut_pangenome.vcf.gz | awk 'length($5) > length($4)' | wc -l```

#### ANSWER: the pangenome contains 12012 biallelic SNPs and 1766 biallelic INDELs, of which 947 are indertions and 827 are deletions.

Of course, these variants needs to be filtered and validated for downstream analyses, but this can give us an idea of the variability among the individuals included in the graph. You can also look at a particular variant with ```SequenceTubeMap``` by chunking the graph around the variant coordinates.

_____

## 5. Mapping of short-reads data with vg giraffe

We will now use the pangenome as a reference for read mapping using the fast short-read mapper ```vg giraffe```<sup>7</sup>.

In the folder ```4_short_read_data/``` you will find forward and reverse fastq files for a single Zebra finch individual [SRR16569049](https://www.ncbi.nlm.nih.gov/sra/SRX12771087[accn]). 
They have been pre-processed as follows:
1. The data have been downloaded from NCBI with ```fasterq-dump```
2. Read statistics were checked with ```fastqc```. No adaptors were found.
3. Read have been aligned to the pangenome with ```vg giraffe``` and surjected with ```vg surject``` (see below for the explanation). Only those mapping to chr22 were extracted.
4. Read have been aligned to the linear reference with ```bwa mem``` (see below for the explanation). Only those mapping to chr22 were extracted.
   
### 5.1 Align the reads with vg giraffe

**RUN:**
```
vg giraffe -t 4 -p \
	--gbz-name 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.gbz \
	-m 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.min  \
    -d 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.dist \
	-f 4_short_reads/SRR16569049_1_chr22.fastq.gz \
    -f 4_short_reads/SRR16569049_2_chr22.fastq.gz \
    -N wildtype10 \
    > 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam
```

It will take around 3.5 minutes to run using 4 threads and 8 GB RAM (see log file below).
Copy and paste the outputs from the [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/5.1_vg_giraffe) folder if there is no time to run the command:
```
cp reference_data/5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam
```

#### The inputs

* ```--gbz-name```: the ```filtered .gbz graph``` generated bu the ```MC pipeline```.
* ```-m``` and ```-d```: the two indexes ```.min``` and ```.dist``` generated by the ```MC pipeline```.
* ```-f```: raw reads to align (forward and reverse)
* ```-N```: the sample name

#### The output

The output of vg giraffe is a ```GAM``` file, which stores alignments in ```vg```’s native binary format.

```vg giraffe``` log (stdout):

```
Using fragment length estimate: 488.569 +/- 128.728
Mapped 1722530 reads across 4 threads in 203.879 seconds with 7.92238 additional single-threaded seconds.
Mapping speed: 2091.87 reads per second per thread
Used 791.928 CPU-seconds (including output).
Achieved 2175.11 reads per CPU-second (including output)
Used 3586241473540 CPU instructions (not including output).
Mapping slowness: 2.08196 M instructions per read at 4528.5 M mapping instructions per inclusive CPU-second
Memory footprint: 0.339954 GB
```

____

#### QUESTION 9: *how many read have aligned?*

To answer this question we need to generate statistics for the ```GAM``` file. We can use the ```vg toolkit``` as follows.

**RUN:**
```vg stats --threads 4 -a 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam > 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam.stats```
**RUN:**
```cat 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam.stats```

```
Total alignments: 1722530
...
Total aligned: 950821
...
```

In the vg stats output, ```total alignments``` are acqually the number of available reads (i.e. possible alignments), while ```total aligned``` are the reads that actually aligned. We have aligned the 55% of the reads, remember we only have the first 1 Mbp of chr22 in the pangenome, while the reads were samples as those aligning to the whole chr22. 

#### ANSWER: we have aligned 950 K paired reads

___

**OPTIONAL - not needed for this excercise**

During the run, ```vg giraffe``` sometimes prints to stdout progress messages like the following:
```Using fragment length estimate: 488.569 +/- 128.728```
These are fragment length estimates that ```vg giraffe``` calculates on the fly (**OUR CASE**).

“In the folder ```5.2_bwa_mem/``` you can find the outputs from the steps below. Even though we won’t need them, it’s important to know how to proceed when ``` vg giraffe``` complains about the insert size.

Sometimes you can get a warning like this one:
```
warning[vg::giraffe]: Encountered 100000 ambiguously-paired reads before finding enough
                      unambiguously-paired reads to learn fragment length distribution. Are you sure
                      your reads are paired and your graph is not a hairball?
warning[vg::giraffe]: Finalizing fragment length distribution before reaching maximum sample size
                      mapped 6 reads single ended with 100000 pairs of reads left unmapped
                      mean: 0, stdev: 1
warning[vg::giraffe]: Cannot cluster reads with a fragment distance smaller than read distance
                      Fragment length distribution: mean=0, stdev=1
                      Fragment distance limit: 2, read distance limit: 200
warning[vg::giraffe]: Falling back on single-end mapping
Using fragment length estimate: 0 +/- 1
```

Giraffe performes single-end mapping by defaults. When aligning paired-end reads it assumes that the first batch of reads are representative of the whole run, but if the input file is sorted, for example, this is not the case (see: https://github.com/vgteam/vg/wiki/Giraffe-best-practices). Therefore, ```vg giraffe``` can’t estimate the mean fragment length and standard deviation for read pairing and falls back to single-end mapping. If this happens, align the reads to a linear reference and calculate those statistics using ```Picard```. 

1. Index the linear reference

```bwa index 1_fasta_files/bTaeGut7_mat_chr22_1Mb_NC_133047.1.fasta```

2. Align the reads to the linear reference with ```Bwa mem```
```
bwa mem -t 2 -M \
-R "@RG\tID:SRR16569049\tPL:ILLUMINA\tSM:wildtype10\tPU:SRR16569049\tLB:SRR16569049" \
1_fasta_files/bTaeGut7_mat_chr22_1Mb_NC_133047.1.fasta \
4_short_reads/SRR16569049_1_chr22.fastq.gz \
4_short_reads/SRR16569049_2_chr22.fastq.gz | \
samtools view -@ 2 -bS |
samtools sort -@ 2 \
-o 5.2_bwa_mem/SRR16569049_bwa_mem_chr22.sort.bam
```

This alingments takes around 5.6 minutes to complete.

Differently to ```vg giraffe```, ```BWA mem``` outputs secondary alignments by default. Remember to remove them for downstream analyses or for comparisons with the ```vg``` alignments.

2. Check the stats for the linear alignment
```samtools flagstats -@ 4 5.2_bwa_mem/SRR16569049_bwa_mem_chr22.sort.bam 1>  5.2_bwa_mem/SRR16569049_bwa_mem_chr22.sort.bam.flagstats.out```

3. Estimate the insert size statistics

```picard CollectInsertSizeMetrics I=5.2_bwa_mem/SRR16569049_bwa_mem_chr22.sort.bam O=5.2_bwa_mem/insert_metrics.txt H=5.2_bwa_mem/insert_hist.pdf M=0.5 VALIDATION_STRINGENCY=SILENT```
```grep -A1 "^MEDIAN_INSERT_SIZE" 5.2_bwa_mem/insert_metrics.txt | tail -n1 | awk -F'\t' '{printf("--fragment-mean %s --fragment-stdev %s\n",$6,$7)}'```

The output is:
```--fragment-mean 392.291361 --fragment-stdev 264.240795```

4. feed the fragment mean size and standard deviation in ```vg giraffe```
```
vg giraffe -t 4 -p \
	--gbz-name 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.gbz \
	-m 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.min  \
    -d 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.dist \
	-f 4_short_reads/SRR16569049_1_chr22.fastq.gz \
    -f 4_short_reads/SRR16569049_2_chr22.fastq.gz \
    -N wildtype10 \
    --fragment-mean 392.291361 --fragment-stdev 264.240795 \
    > 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22_frament_size.gam
```
___

### 5.2 Project graph alignments onto a linear reference with ```vg surject```

```vg giraffe``` produces graph-coordinate alignments (GAM). Most downstream tools such as ```samtools```, ```Picard```, ```bcftools``` and ```GATK```, require linear reference coordinates and therefore linear alignment files (BAM). With surjection you can projects each graph alignment onto a chosen reference path (e.g., the backbone reference) and generate a BAM file for downstreama anlysis. 

#### 5.1 Retrieve the reference paths we will suject the alignments onto

Extract the reference paths (a single "path" in our case) from the path file we generated before. It's fine to use those from the ```clip graph``` since the backbone reference paths are always unclipped.

**RUN:**
```grep "bTaeGut7_mat#0#chr22" 3_stats_and_viz/bTaeGut_pangenome.og.paths > 5.1_vg_giraffe/bTaeGut_pangenome.og.REF.paths```

**RUN:**
```cat 5.1_vg_giraffe/bTaeGut_pangenome.og.REF.paths```

```bTaeGut7_mat#0#chr22```

#### 5.2 run ```vg surject```

**RUN:**
```
vg surject \
	--threads 4 \
	--xg-name 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.xg  \
	-F 5.1_vg_giraffe/bTaeGut_pangenome.og.REF.paths \
    -R "@RG\tID:SRR16569049\tPL:ILLUMINA\tSM:wildtype10\tPU:SRR16569049\tLB:SRR16569049" \
    --interleaved \
    --bam-output \
	5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam \
	> 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam
```

The sujection will take a couple of minutes to complete.
Copy and paste the outputs from the [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/5.1_vg_giraffe) folder if there is no time to run the command:
```
cp reference_data/5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam 
```
#### The inputs
* ```--xg-name```: ```filtered graph``` in ```.xg``` format. We already ave it since we asked MC to generate an ```.xg``` index for the ```clip``` and ```filter``` graph
* ```-F``` : reference paths file
* ```-R``` : read groups information fow downstream analyses. See below.
* ```--interleaved```: to indicate the GAM file contains interleaved paired-end reads
* ```--bam-output```: to indicate we want a bam file as output
* the ```GAM``` input file

For downstream analyses and post-processing, such as  ```picard MarkDuplicates ``` or  ```GATK```, read groups are extremely important to distinguish sets of reads that were generated from a single run of a sequencing instrument, being able to track sequencing library artifacts and technical variations that arise during different sequencing runs. More information of read groups can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).

However, NCBI ofter removes the read group information from the fastq files, so we had to make up unique names.

### The output

The output is a ```BAM``` file referenced to the reference paths we provided (the backbone reference).

#### Sort the bam:

**RUN:**
```samtools sort -@ 4 -o 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam```

#### Reheader the bam:

Remember to rename the ```.bam``` with "standard" chromosome names. The ```PanSN-spec``` sequence naming is not compatible with all the other softwares (e.g. ```mapDamage2```)

**RUN:**
```samtools reheader -c "sed s/bTaeGut7_mat#0#//g" 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam > 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam```

Check before and after.

**Before**

**RUN:**
```samtools view -H 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam```
```
@HD	VN:1.5	SO:coordinate
@SQ	SN:bTaeGut7_mat#0#chr22	LN:1000000	M5:df4418ac75defc621ef63b32e2294cd7
@RG	ID:@RG\tID:SRR16569049\tPL:ILLUMINA\tSM:wildtype10\tPU:SRR16569049\tLB:SRR16569049	SM:wildtype10
@PG	ID:0	PN:vg
@PG	ID:samtools	PN:samtools	PP:0	VN:1.22.1	CL:samtools sort -@ 4 -o 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.22.1	CL:samtools view -H 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam
```

**After**

**RUN:**
```samtools view -H 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam```

```
@HD	VN:1.5	SO:coordinate
@SQ	SN:chr22	LN:1000000	M5:df4418ac75defc621ef63b32e2294cd7
@RG	ID:@RG\tID:SRR16569049\tPL:ILLUMINA\tSM:wildtype10\tPU:SRR16569049\tLB:SRR16569049	SM:wildtype10
@PG	ID:0	PN:vg
@PG	ID:samtools	PN:samtools	PP:0	VN:1.22.1	CL:samtools sort -@ 4 -o 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.bam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.22.1	CL:samtools reheader -c sed s/bTaeGut7_mat#0#//g 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.bam
@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.22.1	CL:samtools view -H 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam
```

This command left only "chr22" as chromosome name!

___

#### QUESTION 10: *how many reads remains after surjection?*

To answer this question we need to generate statistics for the BAM file. We can use ```samtools flagstat``` as follows.

**RUN:**
```samtools flagstats -@ 4 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam 1> 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam.flagstats.out``` 
**RUN:**
```cat 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sort.reheadered.bam.flagstats.out```

```
...
1722530 + 0 in total (QC-passed reads + QC-failed reads)
...
872688 + 0 mapped (50.66% : N/A)
```

#### ANSWER: we have aligned 872 K paired reads

___

#### QUESTION 11: *did we lost some read during surjection?*

Aligned reads in the ```GAM``` file: **950821**
Aligned reads in the ```BAM``` file: **872688**

#### ANSWER: 78133 reads that were mapped in the graphical alignment ``GAM``, are now stored in the ``BAM`` file as unmapped. These reads couldn't be projected one the linear path we chose because they probably aligned exclusively to an alternative path with no embedding on your chosen paths, i.e. true novel insertions), it only flanks/partially overlaps the path so no continuous projection exists, or the alignment is too complex/ambiguous for a consistent projection.

*New computational methods and file formats other than the linear binary alignment map (``BAM``) and variant call format (``VCF``) need to be developed to overcome this limitation and represent all the information embedded in the graph<sup>1</sup>.*

___

#### 5.3 Visualize aligned reads with ```SequenceTubeMap```<sup>6</sup>

```SequenceTubeMap```<sup>6</sup> can be used to visualize a ```GAM``` file. 

As before, let's try to use the [online demo](https://vgteam.github.io/SequenceTubeMap/). If it doesn't work, you can look at my screen and the screenshots and videos I uploaded in the folder 6_vg_giraffe_viz/sequenceTubeMap:

First, we will chunk the ```GAM``` and the ```filter``` graph used for the alignment at the same coordinates we used to chunk the ```clip``` graph before. 

##### 5.3.1 Sort and index ```GAM```

**RUN:**
```vg gamsort -t 4 -i 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam.gai 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.gam  > 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam```

This will take around 3 minutes to run.
Copy and paste the outputs from the [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/5.1_vg_giraffe) folder if there is no time to run the command:
```
cp reference_data/5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam
cp reference_data/5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam.gai 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam.gai
```

##### 5.3.2 Chunk the sorted ```GAM``` and the ```filter``` graph

**RUN:**
```vg chunk -t 4 -c 1 -p bTaeGut7_mat#0#chr22:0-100000 -x 2_bTaeGut_pangenome/bTaeGut_pangenome.d1.xg -a 5.1_vg_giraffe/SRR16569049_vg_giraffe_chr22.sorted.gam -g -O vg --prefix 6_vg_giraffe_viz/bTaeGut_pangenome > 6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.vg```

*Ignore the warning: "warning[vg chunk]: the vg-protobuf format is DEPRECATED. you probably want to use PackedGraph (pg) instead"*

**Inputs:**
* ```-p```: the coordinates to chunk
* ```-x```: the ```.xg``` graph
* ```-a```: the sorted and indexed ```.gam``` file
* ```-g```: tells the program to chunk both the graph and the ```.gam```
* ```-O``` specifies the file name for the chunked graph (```.vg``` format)
* ```--prefix```: specifies the chunked ```.gam``` path and prefix
  
**Outputs:**
* ```6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.gam```: the chunked ```.gam```
* ```6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.vg```: the chunked ```filtered``` graph

##### 5.3.3 Index the chunked ```.gam``` and ```filter``` graph

**RUN:**
```vg convert -t 4 -x 6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.vg > 6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.xg```
**RUN:**
```vg gamsort -t 4 -i 6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam.gai 6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.gam  >  6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam ```

This will take a couple of minutes to run.
Copy and paste the outputs from the [reference_data/](https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/tree/main/reference_data/6_vg_giraffe_viz) folder if there is no time to run the command:
```
cp reference_data/6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam 6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam
cp reference_data/6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam.gai 6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sort.sorted.gai
```

##### 5.3.4 Upload the files in the [online demo](https://vgteam.github.io/SequenceTubeMap/). 

You can find the files in this github repository, download it directly from her the the computer:
* ```6_vg_giraffe_viz/6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.vg```
* ```6_vg_giraffe_viz/6_vg_giraffe_viz/bTaeGut_pangenome.d1.chunk.100k.xg```
* ```6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sorted.gam```
* ```6_vg_giraffe_viz/bTaeGut_pangenome_0_bTaeGut7_mat#0#chr22_0_100016.sort.sorted.gai```

Follow the instructions:
* Go to the [SequenceTubeMap demo page](https://vgteam.github.io/SequenceTubeMap/). 
* Select "Custom" from the "Data" drop down menu > Click on "Configure Tracks" > click the "+" button > leave "graph" but change "mounted" with "upload" > select the .xg file from the Download folder > close using the "x" in the upper right corner
* Click "+" > change "graph" with "read" > change "mounted" with "upload" > select the .gam file from the Download folder
  
If the online demo doesn't work, you can look at my screen and the screenshots and videos I uploaded in the folder ```6_vg_giraffe_viz/sequenceTubeMap/```:





##other things


### 4.5 align reads against the linear reference

bwa index 1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta

bwa mem -t 16 -M \
-R "@RG\tID:SRR16569049\tPL:ILLUMINA\tSM:wildtype10\tPU:SRR16569049\tLB:SRR16569049" \
1_fasta_files/bTaeGut7_mat_chr22_NC_133047.1.fasta \
4_short_reads/SRR16569049_1_chr22.fastq.gz \
4_short_reads/SRR16569049_2_chr22.fastq.gz | \
samtools view -@ 16 -bS |
samtools sort -@ 16 \
-o 5.2_bwa_mem/SRR16569049_bwa_mem_chr22.bam

BWAmem outputs secondary alignments, remember to remove also those to match BWAmem and vg safari

### 4.5 generate statistics for the linear BAM

samtools flagstats -@ 32 5.2_bwa_mem/SRR16569049_bwa_mem_chr22.bam 1> 5.2_bwa_mem/SRR16569049_bwa_mem_chr22.bam.flagstats.out 

#to get insert size

picard CollectInsertSizeMetrics I=5.2_bwa_mem/SRR16569049_bwa_mem.sort.bam O=5.2_bwa_mem/insert_metrics.txt H=5.2_bwa_mem/insert_hist.pdf M=0.5 VALIDATION_STRINGENCY=SILENT
grep -A1 "^MEDIAN_INSERT_SIZE" 5.2_bwa_mem/insert_metrics.txt | tail -n1 | awk -F'\t' '{printf("--fragment-mean %s --fragment-stdev %s\n",$6,$7)}'
--fragment-mean 429.131808 --fragment-stdev 231.084784

### 4.5 remove secondary and supplementary alignments and unmapped reads

While for vg giraffe and surject you need to specify that you also want supplementary and secondary alignments, BWA meme does it by default, so we need to remove those.

samtools view -@ 16 -bF 2308  5.2_bwa_mem/SRR16569049_bwa_mem.sort.bam | samtools sort -@ 16 -o  5.2_bwa_mem/SRR16569049_bwa_mem_mapped.sort.bam
samtools flagstats -@ 32 5.2_bwa_mem/SRR16569049_bwa_mem_mapped.sort.bam 1> 5.2_bwa_mem/SRR16569049_bwa_mem_mapped.sort.bam.flagstats.out 





# References

1. Secomandi, Simona, et al. "Pangenome graphs and their applications in biodiversity genomics." Nature Genetics 57.1 (2025): 13-26.
2. Hickey, Glenn, et al. "Pangenome graph construction from genome alignments with Minigraph-Cactus." Nature biotechnology 42.4 (2024): 663-673.
3. Li, Heng, Xiaowen Feng, and Chong Chu. "The design and construction of reference pangenome graphs with minigraph." Genome biology 21.1 (2020): 265.
4. Armstrong, Joel, et al. "Progressive Cactus is a multiple-genome aligner for the thousand-genome era." Nature 587.7833 (2020): 246-251.
5. Guarracino, Andrea, et al. "ODGI: understanding pangenome graphs." Bioinformatics 38.13 (2022): 3319-3326.
6. Beyer, Wolfgang, et al. "Sequence tube maps: making graph genomes intuitive to commuters." Bioinformatics 35.24 (2019): 5318-5320.
7. Sirén, Jouni, et al. "Pangenomics enables genotyping of known structural variants in 5202 diverse genomes." Science 374.6574 (2021): abg8871.
