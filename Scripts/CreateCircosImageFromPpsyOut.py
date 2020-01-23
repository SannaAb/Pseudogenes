#!/usr/bin/python
import sys
import os
import pysam
import argparse
import pandas as pd
from math import *
from decimal import *
 
# Important it increases 200 nt on each site of the insert site

def parseArgs():
    parser = argparse.ArgumentParser(description='Creates Circos plot from the Ppsy output. Complement to the coverage plots created in GVIZ')    
    parser.add_argument('-S', dest='Sample', help='Sample name, descides the prefix of your output plot for circos, make sure it is unique. It will create a folder with this name and output all the circos outs into this folder', required=True)
    parser.add_argument('-rawBam', dest='baminput',help='Bamfile containing alignment performed by a splice aware aligner, need index in the same folder (required)',required=True)
    parser.add_argument('-ParentgeneCoordinate', dest='parentcoord',help='Coordinates of the parent gene (write as chrN:start-end)',required=True) 
    parser.add_argument('-Insertsite', dest='insertsite',help='Coordinates of the insert of the pseudogene (write as chrN:start-end)',required=True)
    parser.add_argument('-GeneName', dest='genename',help='Gene name of the parent gene, this is used for grepping the coordinates of the exons in the exon database, therefore it is very important that you type it exactly the way that the exon coords are typed',required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def databases():
    """
    If the Ppsy environment is loaded this is the exact path to the exon coord, if you install it in another way you need to change it concordantly 
    """
    exoncoords = os.path.abspath(__file__).split("bin/Ppsy.py")[0] + "HG19_databases/Exon_coord_hg19_refgene.bed"
    #exoncoords = "~/.conda/envs/Ppsy/HG19_databases/Exon_coord_hg19_refgene.bed"
    return(exoncoords)

def CreatesIntermediateFilesCircos(Sample,baminput,parentcoord,insertsite,genename, exoncoords):
    """
    This part creates all intermediate files for plotting with circos
    Creates the karyotype for circos
    Coverage circos format over Parent gene all
    Coverage circos format over Parent gene splice reads 
    Coverage circos format over fusion 
    Create link between the exons
    """
    command = "mkdir -p %s" %Sample 
    os.system(command)
    print "Setting up outputfiles"
    genomefile = Sample+"/"+ Sample + "_" +  genename + ".genome"
    exactcoordbed = Sample+"/"+ Sample + "_" +  genename + ".exoncoords.bed"
    insertcoordbed = Sample+"/"+ Sample + "_" +  genename + ".insertcoord.bed"
    covoverbam_circosformat = Sample+"/"+ Sample + "_"+genename + ".Cov_overBam.depth_circos"
    covoversplitbam_circosformat = Sample+"/"+ Sample + "_"+genename + ".Cov_overSplitBam.depth_circos"
    covoverbamInsert_circosformat = Sample+"/"+ Sample + "_"+genename + ".Cov_overBamInsertsite.depth_circos"
    karyotype = Sample+"/"+ Sample + "_" + genename + "_Karyotype"
    links = Sample+"/"+ Sample + "_LinksbetweenExons.txt"
    print "Creating Karyotype for exons of the gene"
    command = "grep -w %s$ %s > %s" %(genename,exoncoords,exactcoordbed)
    os.system(command)
    with open(karyotype, "w") as karyotypeExons:
            with open(exactcoordbed, "r") as exoncoordread:
                with open(insertcoordbed, "w") as incordfile: 
                    # First print the InsertSite to the karyotype, As default we are adding 200 nt on each side of the fusion to increase the size of the plot 
                    inschrom = insertsite.split(":")[0]
                    instart = int(insertsite.split(":")[-1].split("-")[0])-200
                    instend = int(insertsite.split(":")[-1].split("-")[-1])+200
                    inssecir = insertsite.replace(":","-")
                    print >> karyotypeExons, "chr - %s %s %s %s InsertSite" %(inssecir,inssecir,instart-instart,abs(instend-instart))
                    # Print the same info to a bed file for extraction but here use the regular coordinates but with the additional 200 on each side
                    print >> incordfile, "%s\t%s\t%s\t%s" %(inschrom, instart, instend,insertsite)
                    for line in exoncoordread:
                        line = line.strip()
                        chromosome = str(line.split("\t")[0])
                        startexoncoord = int(line.split("\t")[1])
                        endexoncoord = int(line.split("\t")[2])
                        identifier = line.split("\t")[3]
                        print >> karyotypeExons, "chr - %s %s %s %s exon" %(identifier,identifier,startexoncoord-startexoncoord, abs(endexoncoord-startexoncoord))
    print "Creating genome file for bedtools"
    command = "samtools view %s -H | grep ^@SQ | cut -f 2,3 | sed 's/SN://g' | sed 's/LN://g' > %s" %(baminput,genomefile)
    os.system(command)
    print "Extract from rawbam, coverage over Parent gene"
    bamparentgeneCigarN = Sample+"/"+ Sample + "_" +  genename + ".ParentGeneCigarFilt.bam.tmp"
    samfile = pysam.AlignmentFile(baminput, "rb")
    with pysam.AlignmentFile(bamparentgeneCigarN, "wb", header = samfile.header) as bam_out:
        for read in samfile.fetch(parentcoord.split(":")[0],int(parentcoord.split(":")[-1].split("-")[0]),int(parentcoord.split(":")[-1].split("-")[-1])):
            if "N" in read.cigarstring:
                bam_out.write(read)
    command = "samtools index %s" %bamparentgeneCigarN
    os.system(command)
    # Get coverage of the ordinary bam, this might break due to large amount of data 
    command = "coverageBed -sorted -b %s -a %s -d -split -g %s | awk 'OFS=\" \" {print $4,$5-1,$5-1,$6}' > %s" %(baminput, exactcoordbed,genomefile,covoverbam_circosformat)
    # Turning this off for now! As we have it and it is slow 
    os.system(command) 
    print "Extract from rawbam, coverage over Fusion" 
    # Obs we cannot have : in the identifier so these are replaced with - 
    command = "coverageBed -sorted -b %s -a %s -d -split -g %s | awk 'OFS=\" \" {print $4,$5-1,$5-1,$6}' | sed 's/:/-/g' > %s" %(baminput,insertcoordbed,genomefile,covoverbamInsert_circosformat)
    # Turning this off for now! As we have it and it is slow
    os.system(command)
    print "Extract from SpliceRead bam"
    command = "coverageBed -sorted -b %s -a %s -d -split -g %s | awk 'OFS=\" \" {print $4,$5-1,$5-1,$6}' > %s" %(bamparentgeneCigarN, exactcoordbed,genomefile,covoversplitbam_circosformat)
    os.system(command)
    print "Creating the link file"
    df_karyotype = pd.read_csv(karyotype, sep = " ", names=['chr', '-','identifier', 'identifier_2', 'start', 'end', 'color'])
    exonends = df_karyotype[df_karyotype['color'] == 'exon'][['identifier','end']][:-1] # get only exons from karyotype, then get only the identifier and end column and finally remove the last exons as you dont want a link from this one
    exonstarts = df_karyotype[df_karyotype['color'] == 'exon'][['identifier','start']][1:] # same as above but here we drop the first link because we dont want a start link from the first exon
    exonslinkdf = pd.concat([exonends.reset_index(drop=1).add_suffix('_1'), exonstarts.reset_index(drop=1).add_suffix('_2')], axis=1).dropna(axis=1,how="any")[['identifier_1','end_1', 'end_1','identifier_2','start_2','start_2']]
    exonslinkdf.to_csv(links, sep = " ", header =False, index = False)
    
    print "Removing tmpfiles"
    os.remove(bamparentgeneCigarN)
    bai=bamparentgeneCigarN+".bai"
    os.remove(bai)
    return(karyotype,covoverbam_circosformat,covoversplitbam_circosformat,covoverbamInsert_circosformat,links)


def PlottingCircos(Sample,genename,karyotype,covoverbam_circosformat,covoversplitbam_circosformat,covoverbamInsert_circosformat,links,insertsite):
    """
    Here we plot the circos plot by creating the config file with paths to all the preprocessed files
    """
    # Calculate the total length of the karyotypes to obtain the angle ofset of the plot
    totallengths = 0 
    linenr = 0 
    with open(karyotype, "r") as f: 
        for line in f: 
            if linenr == 0:  # this is the first line that is the insertsite length
                insertsitelength = int(line.split(" ")[5])
            totallengths += int(line.split(" ")[5])
            linenr += 1
    fractionfordegrees = Decimal(360)/Decimal(totallengths)
    angle=90-insertsitelength/2*fractionfordegrees
    # for setting the radious for the inser site
    inssecir = insertsite.replace(":","-")
    ## Create the config file
    configurationfile = Sample+"/"+ Sample + "_" +  genename + ".conf"
    # Image files
    outputpicture = Sample+"/"+ Sample + "_" +  genename + ".png"
    outputpicturesvg = Sample+"/"+ Sample + "_" +  genename + ".svg"
    # Calculate absolute depth 
    depthlist = []
    print "Calulating the max depth for scaling the coverage over exons"
    with open(covoverbam_circosformat, "r") as d: # looping through the regular coverage from the bamfile over the exons. Save the largest value, this will be set as a max used both for the histogram and the heatmap
        for line in d:
            line=line.strip()
            coveragebamdepth = line.split(" ")[3]
            depthlist.append(coveragebamdepth)
    maxcoverage = max(map(int,depthlist))

    print "Calculating the max deph over the insert site"
    depthlistinsert=[]
    with open(covoverbamInsert_circosformat, "r") as d: # 
        for line in d:
            line=line.strip()
            coveragebamdepth = line.split(" ")[3]
            depthlistinsert.append(coveragebamdepth)
    maxcoverageins = max(map(int,depthlistinsert))*2 # i increase it, otherwise it takes over the image

    print "Creating the config file"
    with open(configurationfile, "w") as conf:
        print >> conf, """

# Plot the Circos figure 

# circos.conf

karyotype = %s

<<include colors_fonts_patterns.conf>>
<colors>
exon = 64,64,64
#InsertSite = 204,0,0
InsertSite = 179, 179, 179	

</colors>

chromosomes_radius=%s:1.4r

<ideogram>
<spacing>
default = 0.005r
</spacing>
radius=0.6r
thickness=20p
fill=yes
#label_format = eval(var(label) =~ /^id.*/ ? "" : var(label))
show_label=false #change here?
label_font=default
label_size=30
label_radius=0.6r
label_parallel=yes
</ideogram>

<image>
file*=%s
<<include etc/image.conf>>
#angle_offset* = 90
angle_offset* = %s
</image>

# OBS you are switching the housekeeping conf here. This is due to all error messages
#<<include /jumbo/WorkingDir/B17-006/PseudoScriptDb/housekeeping.conf>> 
<<include etc/housekeeping.conf>>

# RGB/HSV color definition, color lists, location of fonts, fill patterns, included from circos distribution
<<include etc/colors_fonts_patterns.conf>>

show_ticks=yes
show_tick_labels=false
<ticks>
radius=1r
color=black
thickness=5p
<tick>
show_labels=false
spacing=100
size=5p
</tick>
</ticks>

<links>
<link>
file=%s
color=153,0,0
radius=0.95r
bezier_radius=0.7r
thickness=4p
</link>
</links>

<plots>
<plot>
file=%s
type=histogram
r1=0.95r
r0=0.75r
fill_color=black
thickness=2
orientation=in
min=0
max=%s
</plot>

<plot>
file=%s
type=histogram
r1=0.98r
r0=0.6r
fill_color=153,0,0
thickness=2
orientation=in
min=0
max=%s
</plot>

<plot>
file=%s
type=heatmap
r0=1.05r
r1=1.15r
color=ylgnbu-9-seq
min=0
max=%s
</plot>
</plots>

            """ %(karyotype,inssecir,outputpicture,angle,links,covoversplitbam_circosformat,maxcoverage,covoverbamInsert_circosformat,maxcoverageins,covoverbam_circosformat,maxcoverage) 
        
    print "Plotting"
    command = "circos --conf %s" %configurationfile 
    os.system(command)
        
        

def main(Sample,baminput,parentcoord,insertsite,genename):
    """
    Main function, the idea is to take the raw bam, zoom into the positions of interest and create Pseuodogene circosplots
    """
    exoncoords=databases()
    (karyotype,covoverbam_circosformat,covoversplitbam_circosformat,covoverbamInsert_circosformat,links)=CreatesIntermediateFilesCircos(Sample,baminput,parentcoord,insertsite,genename,exoncoords)
    PlottingCircos(Sample,genename,karyotype,covoverbam_circosformat,covoversplitbam_circosformat,covoverbamInsert_circosformat,links,insertsite)
    
if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.Sample,arguments.baminput,arguments.parentcoord, arguments.insertsite, arguments.genename)
