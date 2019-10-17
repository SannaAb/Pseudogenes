# Pseudogenes

## Introduction  
Ppsy finder is pipeline for detecting novel processed pseudogenes using DNA sequencing data (Exomic, Genomic, Targeted gene panels etc). 
Processed pseudogenes are structures that are reintroduced into the genome by retrotransposition. This feature is used by Ppsy finder that detects pseudogene candidates by searching for spliced genes withing the genomic sequencing data. Insert positions of the pseudogene candidates are recorded by linking the pseudogene candidate with softclipping (chimeric reads) and read pair insert sizes (chimeric pairs). The idea of this is dependent on an splice aware aligner that allows softclipping and chimeric read pairs.  

## Installation 

The recommended way of installing the tool is through Anaconda. 

First step is to create a conda evironment for Ppsy which builts upon python 2.7 and load the newly built conda environment using source activate 

```

conda create -n Ppsy python=2.7
source activate Ppsy 

```

When the conda environment for Ppsy is loaded we download the latest version of the tar.ball from the Ppsy github repo and install the downloaded tarball using pip. 

```

wget https://github.com/SannaAb/Pseudogenes/archive/v.0.1.5.tar.gz
pip install v.0.1.5.tar.gz

```
The required python modules (pandas, pysam and psutil) will be installed. 
You can make sure that everything was installed correcly by typing Ppsy.py within your environment.  


### Other dependencies 

Ppsy is also dependent on some tools outside of python. You need to have these tools within your path to make sure that the script works correctly. You can install these tools within your environment using anaconda as well. Make sure that you have source activated your environment first. 

```

conda install -c bioconda bedtools
conda install -c conda-forge -c bioconda samtools bzip2
conda install -c bioconda star
conda install -c bioconda bioconductor-gviz
conda install -c conda-forge ncurses # (?)


````

## How to run 

Ppsy takes either the two paired quality filtered fastq files or the alignment file as input. The parameters are described below 

Tabset {.tabset .tabset-fade .tabset-pills}
## Tab 1
text 1
## Tab 2
text 2
### End tabset



## Structure  
The Ppsy finder consists of two scripts, the Ppsy.py scripts that is the pipeline itself and the script Makeppsyreport.py that outputs a html report for your samples, the results have evidence from both chimeric reads and chimeric pairs. 


## Ppsy Idea 

The main goal of PPsy is to detect inserted processed pseudogenes within human DNA sequencing data. The pipeline is utilizing the fact that processed pseudogenes does not contain any introns. Pseudogene candidates are detected as genes contaning spliced reads across the exon exon junctions. The insert site of the pseuedogene candidates are identified with chimeric read pairs and softclipping. 

The steps of the pipeline is described below 

### Discovering Pseudogene candidates 

The Spliced reads are extracted by screening for the cigarN in the alignment file the resulting reads are saved in an alignment file. Regions with a user defined depth (default 5) is extracted from the newly created alignmentfile, extracted regions with overlapp of atleast 2 exons in an exon coord database are saved as our detected pseudogene candidates. When no pseudogene candidates are found the script is killed. 

### Clipped read extraction 

Clipped reads are reads that are split so one part of the read mapps at one site in the genome and the other part mapps in another site. The Clipped reads are defined in the STAR alignment by the tag SA:. If pseudogenes candidates were detected the clipped reads are extracted and saved in a new alignmentfile. 

### Chimeric read extraction

Chimeric reads are reads with larger insertsize than expected. If pseudogenes candidates were detected the chimeric read pairs are extracted based on a userdefined insertsize treshold, default is 200 000. To small insertsize might increase the amount of false positives. The chimeric read pairs are saved in a new alingment file. 

### Chimeric read overlapp with Pseudogene candidates 

Chimeric reads within the pseudogene candidates are extracted. The extracted chimeric reads are binned together forming anchor evidence. 

* Pseudogene anchor Left start
* Pseudogene anchor Left end
* Fusion Anchor left start
* Fusion Anchor left end

* Pseudogene anchor Right start
* Pseudogene anchor Right end
* Fusion Anchor Right start
* Fusion Anchor Right end

pic? ()

Reads are binned together if they are not further away then a userdefined distance from the starting point. The default distance is 500. Bins are saved if you have enough reads supporting it. The user defines the treshold (default 10, minimum 5). 
The first anchors are called the left anchors, if the following read pairs fusionanchor is within the distance from starting out of the fusion but further away then 5000 bp from the pseuendogene anchor we have hit the right anchors. The right anchors are binned in the same manner as the left ones.

### Clipped reads binning 

The genome coords from the clipped reads are extracted from the clipped alignment file. The clipped reads are saved in ranges if they are within a userdefined distance from the first read (default is 10) for both the first part of the read and the second part of the read. If the amount of user defined supporting reads are enough (default 10, minimum 5) the range is saved. If either part of the read range is within the coordinate as a pseodogenecandidate the softclipping funsionpoint is saved. 

### Combining clipped reads with chimeric read pairs evidence 

The support for the fusions from both chimeric read pairs and clipped reads are combined. The optimal fusion have evidence from both chimeric read pairs (both the left and the right anchor) and the clipped reads. 


### Fusion point annotation 

The detected fusion points (with atleast evidence from the chimeric pairs or the clipped reads) are annotated using annovar. If the fusionpoint is detected using the clipped reads the clipped read start coord sets the annotation of the fusionpoint. If the fusionpoint does not have evidence from the the clipped reads the Fusion left anchor start point is used. 

### Plotting 

The coverage are plotted using GVIZ. 



### Known pseudogene annotation 

The p


## Installation

The easieset way of installing the dependencies are to use conda. 

conda create -n PPsy 

## Dependencies 

## Memory Cons... 

There is no easy way of implementing the total memory consumtion for the script as it is constantly changing to new processes. The easiest way is instead of running the memory profiles which is a python tools that allows you to plot the memory consumtion through time together with the main scripts child processes. 

See examples of how to run it below 

```

mprof run --include-children ../Scripts/Ppsy.py --Method Bam -I SMAD4_chr12_580000_ORD.r1.fqAligned.sortedByCoord.out.bam -S SMAD4_chr12_580000

mprof plot mprofile_20190925095608.dat --title STAR -O STAR_test 

```