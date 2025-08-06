## Pangenome graphs and their application to biodiversity genomics

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

To construct the pangenome we will use the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)<sup>1</sup>. The pipeline can be run as a whole or run step by step. 

#### Main steps of the MC pipeline:

1. A user-selected reference genome is used as the initial backbone
2. The reference is progressively augmented with structural variation from the other genomes by minigraph, a sequence-to-graph aligner, as a graph constructor. The resulting graph is SV only (>50 bp)
3. All assemblies are aligned back to the graph with a minimap2-like algorithm that generates base-level alignments for each reference chromosome separately.
4. A modified version of the reference-free aligner [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) is used to combine the alignments into base-level pangenome graphs that contain variants of all sizes
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
--gfa filter clip full \
--gbz filter \
--giraffe filter \
--vg clip filter full \
--xg clip filter full \
--odgi full \
--viz clip filter full \

```
* ```cactus-pangenome```: calls the MC pipeline 
* ```jobStore```: a directory used by the workflow management system to store intermediate files and job metadata
* ```bTaeGut.seqfile```: the file containing genomes IDs and fasta paths

Multiple flags can be set:

* ```--outDir```: folder where all the outputs will be generated
* ```--outName```: output prefixes
* ```--logFile```: path where to store the log file
* ```--reference```: the reference you want to use as a backbone followed by the other you want to use as reference for the output VCF
* ```--otherContigs``` (optional): tells MC to also consider unassembled scaffolds when generating the graph. They will be grouped in the same subgraph before merging the chromosomes. The name is used-defined. **Here is commented since we only have one chromosome.**
* ```--vcf```: tells MC to generate VCFs files referenced to the ```--vcfReference``` genomes 

In the following tags you can specify on which type of graph you want to operate on. Three type of graphs can be generated from the MC pipeline: 
1. **full graph**: all sequenced are included (default MC output)
2. **clip graph**: sequences not aligned the the SV-only graph produced by Minigraph and other sequences are removed to have a "cleaner" graph. 
3. **filter graph**: used for ```vg giraffe```, contains only nodes trasversed by X number of haplotypes (10% in the HPRC human pangenome)

* ```--filter 1```: removes nodes covered by less that 1 haplotype. We will use 1 since the 10% of our 4 haplotypes would be 0.4 
* ```--gfa clip filter full```: output the graphs in the text-based Graphical Fragment Assembly (GFA) format, which si the default format. This format is typically the most compatible for exchanging graphs between ```vg``` and other pangenome tools
* ```--gbz filter```: generate the .gbz graph file for these type of graphs. The filter grah in .gbz format is needed for ```vg giraffe```.
* ```--giraffe filter```: generate ```vg giraffe``` indexes for the filter graph (default)
* ```--vg clip filter full```:  generate the graphs in the Variation Graph (VG, .vg) format, usefull to use with the vg toolkit
* ```--xg clip filter full```: generate the .xg index file for these type of graphs. It's a compressed, indexed representation of a variation graph, specifically optimized for fast path and graph traversal operations
*  ```--og full```: generate the graph in the ```odgi``` format  ```.og```. Usefull when running operations using ```odgi```
* ```--viz full```: generate an ```odgi viz``` 1D .png for each chromosome. ```odgi``` works better with full graphs, the presence off all sequences doesn't hinder the visualization.

If you have multiple chromosome you can also add ```--chrom-og```  and ```--chrom-vg``` to generate ```.og``` and ```.vg``` files for each chromosome. If you forget to set some of these flags don't worry, you can always generate ```.vg```, ```.og```, and other files later on.

## Emphasis

*This text will be italic*
_This will also be italic_

**This text will be bold**  
__This will also be bold__

_You **can** combine them_

## Lists

### Unordered

* Item 1
* Item 2
* Item 2a
* Item 2b
    * Item 3a
    * Item 3b

### Ordered

1. Item 1
2. Item 2
3. Item 3
    1. Item 3a
    2. Item 3b

## Images

![This is an alt text.](/image/sample.webp "This is a sample image.")

## Links

You may be using [Markdown Live Preview](https://markdownlivepreview.com/).

## Blockquotes

> Markdown is a lightweight markup language with plain-text-formatting syntax, created in 2004 by John Gruber with Aaron Swartz.
>
>> Markdown is often used to format readme files, for writing messages in online discussion forums, and to create rich text using a plain text editor.

## Tables

| Left columns  | Right columns |
| ------------- |:-------------:|
| left foo      | right foo     |
| left bar      | right bar     |
| left baz      | right baz     |

## Blocks of code

```
let message = 'Hello world';
alert(message);
```

## Inline code

This web site is using `markedjs/marked`.

# References

1. Hickey, Glenn, et al. "Pangenome graph construction from genome alignments with Minigraph-Cactus." Nature biotechnology 42.4 (2024): 663-673.


# References

1. Hickey, Glenn, et al. "Pangenome graph construction from genome alignments with Minigraph-Cactus." Nature biotechnology 42.4 (2024): 663-673.
