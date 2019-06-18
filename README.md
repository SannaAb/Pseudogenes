# Pseudogenes

## Introduction  
Ppsy finder is pipeline for detecting novel processed pseudogenes using DNA sequencing data (Exomic, Genomic, Targeted gene panels etc). 
Processed pseudogenes are structures that are reintroduced into the genome by retrotransposition. This feature is used by Ppsy finder that detects pseudogene candidates by searching for spliced genes withing the genomic sequencing data. Insert positions of the pseudogene candidates are recorded by linking the pseudogene candidate with softclipping (chimeric reads) and read pair insert sizes (chimeric pairs). The idea of this is dependent on an splice aware aligner that allows softclipping and chimeric read pairs.  


## Structure  
The Ppsy finder consists of two scripts, the Ppsy.py scripts that is the pipeline itself and the script Makeppsyreport.py that outputs a html report for your samples, the results have evidence from both chimeric reads and chimeric pairs. 

## Installation

The easieset way of installing the dependencies are to use conda. 

conda create -n PPsy 

## Dependencies 

The script is written in python2.7, tools required are:  

*Annovar  
*awk  (?)
*Bedtools  
*Samtools  
*(STAR)  

The script is built within python 2.7, required modules are listed below 

*sys  
*os  
*pandas  
*itertools  
*csv  
*re  
*subprocess  
*glob  
*shutil  
*pysam  
*collections  
*collections import Counter  
*time  
*argparse  
*logging  

Finally the pipeline requires annotation files, the main scripts are based on refgene hg19, the user can change this but in that case make sure that you keep the format. 

*/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Gene_coord_hg19_refgene.bed  
*/jumbo/WorkingDir/B17-006/PseudoScriptDb/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed    
*/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Exon_coord_hg19_refgene.bed  

