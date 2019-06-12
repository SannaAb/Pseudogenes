# Pseudogenes

Scripts for the pseudogene project. 

InsertingTranscriptsandCreatingStructuralVariants.py 
For testing the master scripts, creates structural rearangements based on users request. 


Dependencies:
Annovar
? awk
Bedtools
? Circos 
Samtools 

Python modules (2.7): 
sys
os
pandas
itertools
csv
re
subprocess
glob
shutil
pysam
collections
collections import Counter
time
argparse
logging

database files: 

/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Gene_coord_hg19_refgene.bed
/jumbo/WorkingDir/B17-006/PseudoScriptDb/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed
/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Exon_coord_hg19_refgene.bed
