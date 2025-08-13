## Pangenome graphs and their application to biodiversity genomics

This repository contains the script and files for the pratical part of the course "Pangenome graphs and their application to biodiversity genomics" held the 8th of September 2025 in Ferrara, Italy, as part of the [SIBE summer school](https://sites.google.com/view/sibesummerschool/home-page).

The title of this course recalls our recent pangenomics review (Secomandi et al., 2025)<sup>1</sup> and citations can be find in the following text. 

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

To construct the pangenome we will use the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)<sup>2</sup>. The pipeline can be run as a whole or run step by step. 

#### Main steps of the MC pipeline<sup>1</sup>:

1. A user-selected reference genome is used as the initial backbone
2. The reference is progressively augmented with structural variation from the other genomes by minigraph<sup>3</sup>, a sequence-to-graph aligner, as a graph constructor. The resulting graph is SV only (>50 bp)
3. All assemblies are aligned back to the graph with a minimap2-like algorithm that generates base-level alignments for each reference chromosome separately.
4. A modified version of the reference-free aligner Progressive Cactus<sup>4</sup> is used to combine the alignments into base-level pangenome graphs that contain variants of all sizes
5. Chromosomal graphs are then combined and post-processed to reduce path complexity by collapsing redundant sequences. 

## 1.1 Fasta input files

We will create a pangenome for **chromosome 12** of two different **Zebra finch (*Taeniopygia guttata*)** individuals publicly available on NCBI.

<img src="https://github.com/SimonaSecomandi/SIBE-summer-school-2025-Pangenome-graphs-and-their-applications-in-biodiversity-genomics/blob/main/Zebra%20Finches_Michael_Lawton_Flickr.jpg" alt="drawing" width="300"/> <br/>
*Credits: [Flickr/Michael Lawton](https://www.flickr.com/photos/michaellawton/5712718319)*<br />

The **backbone reference** will be the new T2T reference genome [bTaeGut7.mat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048771995.1/), specifically the maternal haplotype. We will also include the paternal haplotype [bTaeGut7.pat](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048772025.1/) and another chromosome-level individual: [bTaeGut2.hap1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051427915.1/) and [bTaeGut2.hap2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051428105.1/).

Differently from ProgressiveCactus, **MC ignores soft-masked bases**, so there is no need to repeat mask the genomes. For both, hard-masking is not recommended.

## 1.2 bTaeGut.seqfile input file 

The main input file for the MC pipeline is the [bTaeGut.seqfile ]() file. Similarly to Progressive Cactus, this text file contains user defined genome IDs and paths to the corresponding fasta file, one assembly per line. The main difference with the Progressive Cactus imput file is the absence of a guide three (Newick format), which is usually the first row in the Cactus file. 

The file is organized like this:

| GenomeID.hap  | path/to/fasta |
| ------------- |:-------------:|
| bTaeGut7_hap1 | ~/2_fasta_files/bTaeGut7.mat_NC_133037.1_chr12.fasta |
| bTaeGut7_hap2 | ~/2_fasta_files/bTaeGut7.pat_CM109762.1_chr12.fasta |
| bTaeGut2.1 | ~/2_fasta_files/bTaeGut2.hap1_CM121069.1_chr12.fasta |
| bTaeGut2.2 | ~/2_fasta_files/bTaeGut2.hap2_CM121111.1_chr12.fasta |

**IMPORTANT:**
* Divergent haplotypes from the same individuals can be defined with ".1" and ".2" after the genome ID.
* The **backbone reference** (we will use *bTaeGut7.hap1*) can only be a single haplotype and **IT NEEDS TO BE CHROMOSOME-LEVEL**. Its chromosomes will be used as a template to align the other's genomes chromosomes/scaffolds and generate chromosomes graphs that will be then merged. This genome will be the main reference for the output VCF (containing the variants inside the pangenome) and will not appear as sample in the VCF. It will be also used as reference for downsteam analyses. 
* The pipeline currently does not support the presence of both haplotypes of the backbone individual in the form bTaeGut7.1 and bTaeGut7.2. The coordinate space for chromosome compartimentalization can only be based on a single haplotype (the backbone genome). However, the pipeline will work using the IDs **bTaeGut1_hap1** (backbone ref) and **bTaeGut1_hap2** (divergent haplotype of the backbone genome). These will be considered as separate samples, but it's always convenient to include the alternate haplotype of the main reference as you will retain information about their variability for downstream analyses (e.g. read mapping and variant calling). In the output VCF file, bTaeGut1_hap2 will appear as haploid (e.g. 0 instead of 1/0).

* In addition to the backbone reference, one may **specify additional assemblies as reference** when running the MC pipeline. The graph will be referenced to the backbone genome, but the additional genome's paths will serve as a reference for graph decomposition, for example, i.e. the pipeline will generate multiple VCF files, one referenced to the backbone reference and the others to the additional references (see below). The chromosome compartimentalization will still rely on the backbone references's chromosomes.

## 1.3 Running the MC pipeline

Run the following command in the main directory:

```
cactus-pangenome \
./jobStore \
./"$file_name".seqfile \
--outDir bTaeGut_pangenome \
--outName bTaeGut_pangenome \ 
--logFile bTaeGut_pangenome/bTaeGut.log
--reference bTaeGut7_hap1 bTaeGut7_hap2 \
--refContigs chr12 \
#--otherContig chr0others \
--vcf \
--vcfReference bPatFas1_hap1 bPatFas1_hap2 \
--filter 1 \
--gfa clip filter \
--gbz filter \
--giraffe filter \
--vg clip filter \
--xg clip filter \
--odgi clip full \
--viz full \

```
* ```cactus-pangenome```: calls the MC pipeline 
* ```jobStore```: a directory used by the workflow management system to store intermediate files and job metadata
* ```bTaeGut.seqfile```: the file containing genomes IDs and fasta paths

Multiple flags can be set:

* ```--outDir```: folder where all the outputs will be generated
* ```--outName```: output prefixes
* ```--logFile```: path where to store the log file
* ```--reference```: the reference you want to use as a backbone followed by the other you want to use as reference for the output VCF
* ```--refContigs```: the reference chromosomes you want MC to base the chromosomes compartimentalization. In our case it's just ```chr12``` of bTaeGut7_hap1
* ```--otherContigs``` (optional): tells MC to also consider unassembled scaffolds when generating the graph. They will be grouped in the same subgraph before merging the chromosomes. The name is used-defined. **Here is commented since we only have one chromosome.**
* ```--vcf```: tells MC to generate VCFs files referenced to the ```--vcfReference``` genomes 

In the following tags you can specify on which type of graph you want to operate on. Three type of graphs can be generated from the MC pipeline: 
1. **full graph**: all sequenced are included (default MC output)
2. **clip graph**: sequences not aligned the the SV-only graph produced by Minigraph and other sequences are removed to have a "cleaner" graph. 
3. **filter graph**: used for ```vg giraffe```, contains only nodes trasversed by X number of haplotypes (10% in the HPRC human pangenome)

* ```--filter 1```: removes nodes covered by less that 1 haplotype. We will use 1 since the 10% of our 4 haplotypes would be 0.4 
* ```--gfa clip filter```: output the graphs in the text-based Graphical Fragment Assembly (GFA) format, which si the default format. This format is typically the most compatible for exchanging graphs between ```vg``` and other pangenome tools
* ```--gbz filter```: generate the .gbz graph file for these type of graphs. The filter grah in .gbz format is needed for ```vg giraffe```.
* ```--giraffe filter```: generate ```vg giraffe``` indexes for the filter graph (default)
* ```--vg clip filter```:  generate the graphs in the Variation Graph (VG, .vg) format, usefull to use with the vg toolkit
* ```--xg clip filter```: generate the .xg index file for these type of graphs. It's a compressed, indexed representation of a variation graph, specifically optimized for fast path and graph traversal operations
*  ```--og full```: generate the graph in the ```odgi``` format  ```.og```. Usefull when running operations using ```odgi```
* ```--viz full```: generate an ```odgi viz``` 1D .png for each chromosome. ```odgi``` works better with full graphs, the presence off all sequences doesn't hinder the visualization.

If you have multiple chromosome you can also add ```--chrom-og```  and ```--chrom-vg``` to generate ```.og``` and ```.vg``` files for each chromosome. If you forget to set some of these flags don't worry, you can always generate ```.vg```, ```.og```, and other files later on.

## Outputs

## 2. Pangenome evaluation

### 2.1 Statistics 
After the generation of the pangenome, the first thing to do is to check the statistics. This can be done with ```odgi stats``` or ```vg stats```.

```bgzip -d -@ 32 5.1_MC18/5.1_MC18.gfa.gz``` <br />
```odgi build -t 32 -g 5.1_MC18/5.1_MC18.gfa -o 5.1_MC18/5.1_MC18.gfa.og``` <br />
```printf "Generating stats for 5.1_MC18/5.1_MC18.gfa.og"``` <br />
```odgi stats -t 32 -S -i 5.1_MC18/5.1_MC18.gfa.og > 5.1_MC18/5.1_MC18.gfa.og.stats```

**The output:**

You can ```cat``` the output from ```odgi stats``` and look at the content:

```cat 5.1_MC18/5.1_MC18.gfa.og.stats```

| length|nodes|edges|paths|steps|
| ----- |:----:|:----:|:----:|:----:|
| 1661172842|107376702|146402494|0|0

The graph has:
- **Length of X bp:** this is the pangenome graph sequences length, the sum of the lengths of all nodes in the graph.
- X nodes
- X edges 

*Is the pangenome graph bigger than the original reference sequence?* 

**YES!!**

Original size: X vs. Pangenome size: X

*The size of a pangenome graph depends on the genome size of the respective species but is bound to be larger, as it incorporates accessory sequences from other individuals, and it is also influenced by the number and diversity of the individuals contributing to the pangenome as well as by the construction pipeline<sup>1</sup> *

*Of how much the reference was augmented by the other sequences?* 

The other sequences augmented the graph by X bp (X %). 

### 2.2 Subsampling and visualization 

Visualization is important to get an idea of the structure of the graph and the variability among the different genomes.

#### odgi viz<sup>5</sup>

In the MC command we specified to generate ```odgi viz``` graphs and these can be found in the folder XX. Odgi viz was generated on the full graph as explained before and a .png file for each chromosome (one in our particular case) was generated. The input for odgi viz was the -og chromosome object.

Let's look at it. Each line represents a different chromosome' path with their genome.ID in the right side. Each path is coloured when passing through a node and edges representing variantion are represented by black lines in the bottom of the figure.

#### SequenceTubeMap<sup>6</sup>

Another useful tool for visualizing pangenome graphs is SequenceTubeMap. It visualizes the graph in ```.vg``` format using the same linear visualization as ```odgi viz```, but variability among genomes is visualized differently and it can be inspected interactively.

To visualize a specific vg file without uploading it on the webpage, it is possible to launch a server which provides the data to SequenceTubeMap<sup>6</sup>. In this course we will use the web interface for simplicity. 

First, we will chunk the graph in a smaller piece to be able to upload it online (the limit is X Mb).

1. Run this command on the graph: 

```vg chunk -t 64 -c 20 -x 5.1_MC18.d4.xg -p bPatFas1_hap1#0#chr1:10000000-100010000 -O vg > Chr1_10k.d4.vg```

The output would be: ...

2. Prepare the graph for SequenceTubeMap. These commands are included in the prepare_vg.sh script in [SequenceTubeMap repository](https://github.com/vgteam/sequenceTubeMap/tree/master/scripts)

```vg convert "${1}" -x >"${1}.xg```
```vg gbwt -x "${1}" -v "${1%.vg}.vcf.gz" -o "${1}.gbwt``` #do i need this??

You can find it in the folder on this github repository, download it directly from her the the computer.

2. Go to the [sequenceTubeMap demo page](https://vgteam.github.io/sequenceTubeMap/). Select "Custom" from the "Data" drop down menu > Click on "Configure Tracks" > click the "+" button > leave "graph" but change the "mounted" with "upload" > select the file from the Download folder. 

Inspect the graph:
- How many variants do you see?
- ....

## 3. Pangenome-embedded variants

*We looked at the variants inside the pangenome, but how can I look at them in a canonical way and use them for downstream analysis?*

The MC pipeline produces VCF files referenced to the backbone reference and other genomes you specified in the command. You can find information about the VCF format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

The process through which the variants are defined is called **graph decomposition**, the process of breaking down a pangenome graph into smaller, more manageable subgraphs or components (snarls or bubbles).

The file names are: X.vcf and X.vcf

You might have noticed "raw" VCF files. These are those directly outputted by vg deconstruct inside that are then normalized and postprocessed automatically.

Let's look at the VCF referenced to out backbone reference. Run the following command:

bcftools view -H 5.1_MC18.raw.vcf.gz | head


