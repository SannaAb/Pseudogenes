#!/usr/bin/python 

import sys
import os 
import pandas as pd
import itertools
import csv
import re
import subprocess
import glob 
import shutil 
import pysam
import collections
from collections import Counter
import time
import argparse
import logging


# Idea: 
# Extract Split reads CigarNs, Intersect with genes? Pseudogene Candidates!
# Extract Softclipped reads, Tresholds
# Remapp Softclipped Parts of the reads 

# Debugging 
# - Break when there are no evidence after anchor treshold
# - Clean loop that can be called everywhere without break, just moves the files if there is a file. If there is not it does not break. 
# - Created a genome file for creating a seperate chrom info

# ... break because of bedtools memory depletion, should be sorted.  
# Set the chrominfo! Obs! Get it from the header of the alignment file always. 


# 6/2-19 problem with coverage over split bam and bam when creating the circos, have a look at $ where it must end with nothing except the word. You might need to look at the database if this is a problem!  

# 15/2-19 You have a problem with calculating the max for a pseudogne that has the same name as something similair we detected before. This is the reg finding that works FNTB but i also find Churc1-FNTB here it crashes. 


# 9/05 - Time to put in a proper breakpoint of the pseuogene. Skip the breakpoint within the gene itself and focus on the insert site. Use a range from the left and the right anchors to get the breakpoint. Not sure how to handle the pseudogenes coords but maybe they are not that important? Only the structure might be of interest 



# 10/05 - The problem with adding the start and the end coords from reads in a range is that the end of the read might be very far away due to an indel? What do to then?  

# New: 
# Debugging 
# Remove Clipped_Chimeric_Evidence.txt because you already have in in the annotation from annovar 

# Adding new databases for the refgene, i am still running the pseudogenes for the ensemble annotation though. Maybe we can change this? Problem is that the refgene does not contains the info about whether or not it is a pseudogene
 
#  12/6 obs! The next anchor for the chim pairs might be 5000 nucl away within the same gene

# 17/6 this part plots using GVIZ instead 

def parseArgs():
    parser = argparse.ArgumentParser(description='Detects processed pseudogenes by looking at DNA data that splits across the splice junctions')
    parser.add_argument('-I', dest='baminput',help='Bamfile containing alignment performed by a splice aware aligner, need index in the same folder (required)',required=True)
    parser.add_argument('-S', dest='Sample', help='Sample name, descides the prefix of your outputs together with the output folder (required)', required=True)
    parser.add_argument('--pseudoCandidateDepth', dest='Psdepth', help='The minimum depth that supports the splice junctions, these are the resulting processed pseudogene candidates (default 5)', default=5 ,type=int) 
    parser.add_argument('--InsertDistance', dest='insdistance', help='What is the distance from the parent gene where we can have an pseudogene, low distance might increase the amount of detected pseudogenes but will also increase the amount of false positives. A low distance and you might hit inserted pseudogenes in the parent gene itself which is not very likely (default 200 000)', default=200000,type=int)
    parser.add_argument('--ChimericPairDepthTreshold', dest='ChimPairDepthTresh', help='The minimum amount of reads to suppport the chimeric pairs in the left anchor, (min 5), (default 10)', default=10,type=int)
    parser.add_argument('--ChimericPairBinningTreshold', dest='ChimPairBinningTresh', help='When the chimeric pairs are binned into the anchors the binning distance defines the distance for the read to belong to same bin, (default 500)', default=500,type=int)
    parser.add_argument('--ChimericReadDepthTreshold', dest='ChimReadDepthTresh', help='The minimum amount of reads to support the chimeric reads in the fusion site, (min 5) (default 10)', default=10,type=int)     
    parser.add_argument('--ChimericReadBinningTreshold', dest='ChimReadBinningTresh', help='When the chimeric reads are binned into the anchors the binning distance defines the distance for the read to belong to same bin', default=10,type=int)
    parser.add_argument('--MergeChimReadWithChimpairTresh', dest='chimreadpairdistance', help='When we are combining the results from the chimeric pairs and the chimeric reads we combine them if the chimeric read are withing the chimeric pair anchors or the chimeric read are in a user defined distance from the left anchor (default 100)', default=100,type=int)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments


def database(baminput,Sample): 
    #baminput = sys.argv[1] 
    #Sample = sys.argv[2]
    chrominfo = Sample + ".ChromInfo_SortingOrder.txt"
    logging.info('%s\tInput Alignment %s', time.ctime().split(" ")[-2],baminput)
    logging.info('%s\tReading the database files', time.ctime().split(" ")[-2])
    genecoords="/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Gene_coord_hg19_refgene.bed"
    pseudogenecoords="/jumbo/WorkingDir/B17-006/PseudoScriptDb/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed"
    exoncoords="/jumbo/WorkingDir/B17-006/PseudoScriptDb/VersionAnnotationrefgene/Exon_coord_hg19_refgene.bed" # This one wont contain overlapping Coords!
    # Creating chrominfo database for your Alignment file
    command = "samtools view %s -H" % baminput 
    a = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    with open(chrominfo, "w") as chromout: 
        for line in a.stdout:
            line = line.strip()
            if line.split("\t")[0] == "@SQ":
                chrom  = line.split("\t")[1].split(":")[-1]
                size = line.split("\t")[2].split(":")[-1]
                print >> chromout, chrom +"\t"+ str(size)
    # Annovar Databases
    anndb = "/apps/bio/apps/annovar/20150322/humandb"
    annovarscript = "/apps/bio/apps/annovar/20150322/annotate_variation.pl"
    # Lists for cleaning and moving the files    
    Cleaninglist=[]
    Cleaninglist.append(chrominfo)
    # Append to cleaning the chromosome info 
    MovingList = [] 
    Outputfolder=Sample + "_PPsyOut"
    command = "mkdir -p %s" % Outputfolder 
    os.system(command)
    return (chrominfo,genecoords,pseudogenecoords,exoncoords,anndb,annovarscript,Cleaninglist, MovingList, Outputfolder) 
    
def extractingClipped(Sample,baminput, Cleaninglist,MovingList):
    '''
    Softclipped reads are extracted and saved to a new bamfile
    '''
    logging.info('%s\tExtracting Softclipped reads',time.ctime().split(" ")[-2])
    bamClipped = Sample + ".Clipped.bam"
    commandExtractingclipped = "samtools view -h %s | awk '$17 ~ /SA/ || $1 ~ /^@/' | samtools view -bS - > %s" % (baminput, bamClipped)
    os.system(commandExtractingclipped)
    command = "samtools index %s" % bamClipped
    os.system(command)
    MovingList.append(bamClipped)
    MovingList.append(bamClipped+".bai")
    return bamClipped

def extractingchimericReads(Sample, baminput, Cleaninglist, MovingList,insdistance): 
    '''
    Extracting Reads with large insertsite Either if the mate maps in another chromosome or the insertsite is larger than a user defined treshold, the default is 200 000   
    '''
    logging.info('%s\tExtracting chimeric reads, ins > %s or seperate chr mapping', time.ctime().split(" ")[-2], insdistance)

    bamChimericPairs = Sample + ".ChimericReadPairs.bam"
    MovingList.append(bamChimericPairs)
    MovingList.append(bamChimericPairs + ".bai")
    command = "samtools view %s -h | awk 'function abs(x){return ((x < 0.0) ? -x : x)} $7 != \"=\" || abs($4-$8) > %s {print $0}' | samtools view -bS - > %s" %(baminput, insdistance,bamChimericPairs)  # Obs i increased it to 200 000 as some of the candates are longer than 100 000 bp 
    os.system(command)
    command = "samtools index %s" % bamChimericPairs 
    os.system(command)
    return bamChimericPairs 

def ChimericReadsOverlapwithPseudogeneCandidates(Sample, bamChimericPairs,Pseudogenecandidatesbed, Cleaninglist, MovingList,ChimPairDepthTresh,ChimPairBinningTresh): 
    '''
    Here we look for intersect between pseudogene candidates and chimericreads overlap 
    This will add evidence to the pseudogene candidates by looking at overlap between these coordinates and the softclip coordinates 
    '''
    logging.info('%s\tLooking for an overlap between chimeric reads and Pseudogenes candidates',time.ctime().split(" ")[-2])
    PseudogeneCandidateChimbed = Sample + ".ChimericPairs_AmountOfReads.bed"
    Cleaninglist.append(PseudogeneCandidateChimbed)
    samchimpericPairs = pysam.AlignmentFile(bamChimericPairs, 'rb')
    with open(PseudogeneCandidateChimbed, "w") as out: 
        with open(Pseudogenecandidatesbed, "r") as pseudogenecandidatesbedf: 
            # For each line in the bedfile Create a new bedfile for the chimeric pseudogenes fusion sites 
            for line in pseudogenecandidatesbedf:
                line = line.strip()
                chrom = line.split("\t")[0]
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                pseudogenecandidate = line.split("\t")[3]
                PseudogeneCandidateChimbedtmp =Sample+"_" + pseudogenecandidate + ".ChimericPairs_AmountOfReads.bed.tmp"
                Cleaninglist.append(PseudogeneCandidateChimbedtmp)
                with open(PseudogeneCandidateChimbedtmp, "w") as tmpout: 
                    for read in samchimpericPairs.fetch(chrom, int(start), int(end)): 
                        if read.tostring(samchimpericPairs).split("\t")[6] == "=":
                            matechrom = chrom # The fusion is on the same chrom as the pseudogene 
                        else: 
                            matechrom = read.tostring(samchimpericPairs).split("\t")[6] # The fusion is in a sperate chrom
                        try:
                            matepositonslist = samchimpericPairs.mate(read).positions # This part might be super slow, maybe you need to read sort this? We extract the positions of the mate togehter with the positions of the reads  
                            startMate = matepositonslist[0]
                            endMate = matepositonslist[-1]
                            startInGene = read.positions[0]
                            endInGene = read.positions[-1]
                            print >> tmpout, chrom + "\t" + str(startInGene) + "\t" + str(endInGene)+"\t"+matechrom + "\t" + str(startMate) + "\t" + str(endMate) 
                        except ValueError:
                            continue
                # Make sure that the output file is not empty,if it is continue with the next coord 
                if not os.stat(PseudogeneCandidateChimbedtmp).st_size==0:
                    logging.info('%s\tOverlapping coord detected',time.ctime().split(" ")[-2])
                    logging.info('%s\tStart binning the left and right anchors, atleast %s i depth',time.ctime().split(" ")[-2],ChimPairDepthTresh)
                    # Binning using pandas df
                    data = pd.read_csv(PseudogeneCandidateChimbedtmp, sep = "\t", header=None)
                    sorted = data.sort_values([0, 3,1,4]).reset_index(drop=True)
                    if len(sorted) == 1: # If you only have 1 hit the loop will break, 
                        continue 
                    sorted[["diff1", "diff2"]] = sorted.groupby(3)[1,4].diff()
                    #print sorted
                    fusionchroms = []
                    groupedsorted = sorted.groupby(3)
                    for name, group in groupedsorted: 
                        group = group.sort_values([3,4,1]).reset_index(drop=True)
                        PseudogeneanchorLeftstarts = []
                        PseudogeneanchorLeftends= []
                        FusionAnchorleftstarts = []
                        FusionAnchorleftends = []
                        PseudoanchorRightstarts = []
                        PseudoanchorRightends = []
                        FusionAnchorRightstarts = [] 
                        FusionAnchorRightends= [] 
                        if len(group) > 4: # We wont use groups with less than 5 items, this is not enough for the binning anyways, if only one hit it will crash
                            groupCounter = 0 # This is the counter that makes sure that the loop is not in the end, if we have a hit but the loop ends we will loose that hit 
                            for i, row in group.iterrows():
                                groupCounter +=1 
                                if groupCounter != len(group): # If we are not in the final hit 
                                    if not name in fusionchroms:# This is the first hit in the group which set up the lists
                                        PseudogeneanchorLeftstarts.append(row.loc[1])
                                        PseudogeneanchorLeftends.append(row.loc[2])
                                        FusionAnchorleftstarts.append(row.loc[4])
                                        FusionAnchorleftends.append(row.loc[5])
                                        fusionchroms.append(name)
                                    else:
                                        if abs(PseudogeneanchorLeftstarts[0] - int(row.loc[1])) < ChimPairBinningTresh and abs(FusionAnchorleftstarts[0] - int(row.loc[4])) < ChimPairBinningTresh: # if the difference from the previous and the new coords is smaller than the user defined trehsolds for both start sites we save them to the binnning, 500 is default 
                                            PseudogeneanchorLeftstarts.append(row.loc[1])
                                            PseudogeneanchorLeftends.append(row.loc[2])
                                            FusionAnchorleftstarts.append(row.loc[4])
                                            FusionAnchorleftends.append(row.loc[5])
                                        elif abs(PseudogeneanchorLeftstarts[0] - row.loc[1]) > 5000 and abs(FusionAnchorleftstarts[0] - int(row.loc[4])) < ChimPairBinningTresh: # If the range within the pseudogene is larger than 5000 you might have found the second anchor. We are further away for the pseudogene anchor then what we would expect, therefore we are have the right anchor if we are in the same range for the fusion still
                                            PseudoanchorRightstarts.append(row.loc[1])
                                            PseudoanchorRightends.append(row.loc[2])
                                            FusionAnchorRightstarts.append(row.loc[4]) # Fusion coord on the right side is the end coord of the read in the pair outside of the gene. Always think you want to have as much space as possible 
                                            FusionAnchorRightends.append(row.loc[5])
                                        else: 
                                            # Now we cannot continue as we are not within the range of anything. If we have enough for the range we count them. Otherwise we skip them and continue the loop                           
                                            if len(FusionAnchorleftstarts) > ChimPairDepthTresh: 
                                                PseudogeneAnchorleftStart=PseudogeneanchorLeftstarts[0]
                                                PseudogeneAnchorleftend= PseudogeneanchorLeftends[-1]# important, the value is from the previous row. The current row is the faulty, when getting the previus value you need to use the group index 
                                                FusionAnchorleftStart = FusionAnchorleftstarts[0]
                                                FusionAnchorleftEnd = FusionAnchorleftends[-1]
                                                leftstring = str(row.loc[3]) + "\t" + str(FusionAnchorleftStart) + "\t" + str(FusionAnchorleftEnd) + "\t"+str(len(FusionAnchorleftstarts))
                                                if len(PseudoanchorRightstarts) > ChimPairDepthTresh: 
                                                    # This part tells us we have the right anchor as well 
                                                    PseudogeneAnchorrightStart=PseudoanchorRightstarts[0]
                                                    PseudogeneAnchorrightEnd=PseudoanchorRightends[-1]
                                                    FusionAnchorrightStart=FusionAnchorRightstarts[0]
                                                    FusionAnchorrightEnd=FusionAnchorRightends[-1]
                                                    print >> out, pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\t"+str(row.loc[3])+ "\t" + str(FusionAnchorrightStart)+"\t" + str(FusionAnchorrightEnd) + "\t" + str(len(PseudoanchorRightstarts))

                                                else: 
                                                    # We did not detect the second anchor
                                                    print >> out, pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\tNA\tNA\tNA\tNA"
                                                    # We are done with all the 
                                            # We did not find anything for this loop, or we did but we reached the end. Therefore we need to reset everything
                                           
                                            # Left sites are the ones we find first 
                                            PseudogeneanchorLeftstarts=[row.loc[1]]
                                            PseudogeneanchorLeftends=[row.loc[2]]
                                            FusionAnchorleftstarts=[row.loc[4]]
                                            FusionAnchorleftends=[row.loc[5]]
                                            # Then we go to the right side, reset. 
                                            PseudoanchorRightstarts=[]
                                            PseudoanchorRightends=[]
                                            FusionAnchorRightstarts=[]
                                            FusionAnchorRightends=[]                                    
                                else: # We are in the end of the loop
                                    if len(FusionAnchorleftstarts) > ChimPairDepthTresh:
                                        PseudogeneAnchorleftStart=PseudogeneanchorLeftstarts[0]
                                        PseudogeneAnchorleftend= PseudogeneanchorLeftends[-1]# important, the value is from the previous row. The current row is the faulty, when getting the previus value you need to use the group index
                                        FusionAnchorleftStart = FusionAnchorleftstarts[0]
                                        FusionAnchorleftEnd = FusionAnchorleftends[-1]
                                        leftstring = str(row.loc[3]) + "\t" + str(FusionAnchorleftStart) + "\t" + str(FusionAnchorleftEnd\
) + "\t"+str(len(FusionAnchorleftstarts))
                                        if len(PseudoanchorRightstarts) > ChimPairDepthTresh:
                                            # This part tells us we have the right anchor as well 
                                            PseudogeneAnchorrightStart=PseudoanchorRightstarts[0]
                                            PseudogeneAnchorrightEnd=PseudoanchorRightends[-1]
                                            FusionAnchorrightStart=FusionAnchorRightstarts[0]
                                            FusionAnchorrightEnd=FusionAnchorRightends[-1]
                                            print >> out, pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\t"+str(row.loc[3])+ "\t" + str(FusionAnchorrightStart)+"\t" + str(FusionAnchorrightEnd) + "\t" + str(len(PseudoanchorRightstarts))
                                        
                                        else: 
                                            # We did not detect the second anchor
                                            print >> out, pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\tNA\tNA\tNA\tNA"
                                    
    samchimpericPairs.close()        # Close the pysam 
    return PseudogeneCandidateChimbed  #  return the name of the output where we have our fusion coords

def ClippedBamtoCoord(bamClipped,Cleaninglist, MovingList):
    '''
    Converting the chimeric ranges to coordinates
    '''
    logging.info('%s\tExtracting the Clipped read coordinates',time.ctime().split(" ")[-2])
    clippedtmp = bamClipped.split(".bam")[0] + ".sam"
    Cleaninglist.append(clippedtmp)
    clippedbed = bamClipped.split(".bam")[0] + ".bed"
    Cleaninglist.append(clippedbed)
    command = "samtools view %s > %s" %(bamClipped, clippedtmp)
    os.system(command)
    listofreads = [] # Making sure that if you extracted the coords from one read that you dont extract them again!
    with open(clippedtmp, "r") as chim:
        with open(clippedbed, "w") as out:
            for line in chim:
                line = line.strip()
                read = line.split("\t")[0]
                if read not in listofreads:
                    chrom = line.split("\t")[2]
                    chromchim = line.split("\t")[16].split(":")[2].split(",")[0]
                    start = line.split("\t")[3]
                    startchim = line.split("\t")[16].split(":")[2].split(",")[1]
                    print >> out, chrom+":"+start+"\t"+chromchim+":"+startchim
                    listofreads.append(read)               
    return clippedbed 

def clippedreadbinning(clippedbed, Cleaninglist, MovingList,ChimReadDepthTresh,ChimReadBinningTresh): 
    '''
    Making sure that we have enought depth  over the chimeric reads, userdefined treshold. Bin the first mapping coord using bedtools merge, then loop through the hits in the chimeric alignment to add them to the first coord 
    '''
    clippedbedrange = clippedbed.split(".bed")[0] + ".clippedrange.bed"
    Cleaninglist.append(clippedbedrange)
    clippedmapping = clippedbed.split(".bed")[0] + ".clipped.bed.withcov.txt"
    Cleaninglist.append(clippedmapping)
    logging.info('%s\tBinning the clipped reads coords, depth %s',time.ctime().split(" ")[-2], ChimReadDepthTresh)
    command = "cut -f 1 %s | sed 's/:/\t/g' | sort -k1,1 -k2,2n | awk 'OFS=\"\t\" {print $1,$2,$2}' | mergeBed -c 1 -o count -d %s| awk '$4 > 4 {print $1\"\t\"$2\"\t\"$3}' > %s" %(clippedbed, ChimReadBinningTresh,clippedbedrange) # obs! Overlap is counted if they overlap a userdefined treshold away in nucleotides
    os.system(command)
    # Need to now how many lines you have so you start the counting before the for loop breaks below 
    with open(clippedbed) as f:
        for nrLines, l in enumerate(f):
            pass
    with open(clippedmapping, "w") as clippedout: 
        with open(clippedbedrange, "r") as clippedrange: 
            counter = 0 
            for line in clippedrange:
                line = line.strip()
                chrom = line.split("\t")[0]
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                with open(clippedbed, "r") as chimeras: 
                    counter += 1
                    d = {}
                    rownr = 0
                    for l in chimeras:
                        rownr += 1 
                        l = l.strip()
                        chromfirst = l.split("\t")[0].split(":")[0]
                        startfirst = l.split("\t")[0].split(":")[1]
                        chromchim = l.split("\t")[1].split(":")[0]
                        startchim = l.split("\t")[1].split(":")[1]
                        if chromfirst == chrom and end > startfirst > start: 
                            if chromchim in d.keys(): # if i have the chromosome for that range append
                                d[chromchim].append(startchim)
                            else:  # Else create as new chromosome
                                d[chromchim] = [startchim]
                                
                    for key, value in d.iteritems():
                        firstitem = 0 
                        amount = 1 
                        if len(value) > ChimReadDepthTresh: # If you have atleast 10 in depth by default 
                            sortedlist = sorted(map(int,value))
                            for i in range(1,len(sortedlist)): # Remeber the indexes
                                if sortedlist[i] - sortedlist[i-1] < ChimReadBinningTresh: # If the previus item is atleast 20 away from the new one by default
                                    amount += 1 
                                    if i+1 == len(value) and amount > ChimReadDepthTresh: # If you are reaching the end of the list and you have atleast 5 
                                        print >> clippedout, line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i])
                                else: 
                                    if amount > ChimReadDepthTresh: # not in the end of the list but you have 10 by default
                                        print >> clippedout,  line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i-1])
                                    firstitem = i 
                                    amount = 1 
    return clippedmapping
                                                                
def intersectClippedPseudogeneCandidates(Sample, clippedmapping, Pseudogenecandidatesbed, Cleaninglist, MovingList):
    clippedwithpseudogeneoverlap = Sample + ".Clipped_overlapPseudogenes.txt"
    Cleaninglist.append(clippedwithpseudogeneoverlap)
    logging.info('%s\tIntersect the clipped binned read coords with pseudogene candidates',time.ctime().split(" ")[-2])
    with open(clippedwithpseudogeneoverlap, "w") as out: 
        with open(Pseudogenecandidatesbed, "r") as pseudocandidates: 
            for line in pseudocandidates:
                line = line.strip() 
                pseudochrom = line.split("\t")[0]
                pseudostart = int(line.split("\t")[1])
                pseudoend = int(line.split("\t")[2])
                pseudogene = line.split("\t")[3]
                with open(clippedmapping, "r") as clipped: 
                    for cline in clipped: 
                        cline = cline.strip()
                        clipfirstchrom = cline.split("\t")[0] 
                        clipsecondchrom = cline.split("\t")[4] 
                        clipfirststart = int(cline.split("\t")[1] )
                        clipsecondstart = int(cline.split("\t")[5] )
                        clipdepth = cline.split("\t")[3] 
                        clipfirstend = int(cline.split("\t")[2] )
                        clipsecondend = int(cline.split("\t")[6])
                        if pseudochrom == clipfirstchrom: # If first coordinates in the clip overlaps the pseudogene candidate 
                            if pseudostart < clipfirststart < pseudoend or pseudostart < clipfirstend < pseudoend: 
                                print >> out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(pseudogene, pseudochrom, pseudostart, pseudoend, clipsecondchrom, clipsecondstart, clipsecondend, clipdepth)                       
                        elif pseudochrom == clipsecondchrom:       
                            if pseudostart < clipsecondstart < pseudoend or pseudostart < clipsecondend < pseudoend:
                                print >> out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(pseudogene, pseudochrom, pseudostart, pseudoend, clipfirstchrom, clipfirststart, clipfirstend,clipdepth)
    
    if not os.stat(clippedwithpseudogeneoverlap).st_size==0: # The clipped reads with the pseudogenes might be empty and this will break the loop. Therefore we are making sure that the file is not empty before the analysis
        pandasdf=pd.read_csv(clippedwithpseudogeneoverlap, sep = "\t", header=None) # If we have the same fusion in different positions in the pseudogene we merge them 
        pandasdf2= pandasdf.groupby([0,1,2,3,4,5,6], as_index=False)[7].sum()
        pandasdf2.to_csv(clippedwithpseudogeneoverlap,sep="\t", header = False, index=False)
    else: 
        logging.info('%s\tObs the Psedugene clipp overlap is emtpy',time.ctime().split(" ")[-2])
    return clippedwithpseudogeneoverlap

def CombiningClippWithChimericReads(Sample, clippedwithpseudogeneoverlap, PseudogeneCandidateChimbed, Cleaninglist, MovingList,chimreadpairdistance):
    '''
    Adding the evidence to the clipped reads from the chimeric reads 
    '''
    clippedlist=[]
    chimlist=[]
    logging.info('%s\tMerging the overlapping coords from chimeric read pairs and chimeric clipped reads',time.ctime().split(" ")[-2])
    clipandchimevidence = Sample + ".Clipped_Chimeric_Evidence.txt"
    DetectedPseudogenesList = [] # List of Pseudogenes that is evident in both chim reads and pairs. We will use thes pseudogenes for the plots  
    Cleaninglist.append(clipandchimevidence)
    with open(clipandchimevidence, "w") as out: 
        print >> out, "Pseudogene\tP_chrom\tP_start\tP_end\tFus_chrom_left\tFus_start_left\tFus_end_left\tnChimPairs_left\tFus_chrom_right\tFus_start_right\tFus_end_right\tnChimPairs_right\tChimRead_chrom\tChimRead_Start\tChimRead_End\tnChimReads"         
        with open(clippedwithpseudogeneoverlap, "r") as clippedCoord: 
            for cline in clippedCoord:
                cline = cline.strip() 
                pseudogene = cline.split("\t")[0]
                pseudogenechrom = cline.split("\t")[1] 
                pseudogenestart = cline.split("\t")[2]
                pseudogeneend = cline.split("\t")[3]
                cfusionchrom= cline.split("\t")[4]
                cfusionstart = int(cline.split("\t")[5])  
                cfusionend = int(cline.split("\t")[6])   
                cfusiondepth = str(cline.split("\t")[7])
                with open(PseudogeneCandidateChimbed, "r") as chimcoord: 
                    for chimline in chimcoord:
                        chimline = chimline.strip()
                        chimpseudogene = chimline.split("\t")[0]
                        chimpsudochrom = chimline.split("\t")[1]
                        chimpseudostart = int(chimline.split("\t")[2])  # Increase the range with 100 on both sides 
                        chimpseudoend = int(chimline.split("\t")[3])  
                        chimfuschrom = chimline.split("\t")[4]
                        chimfusstartleft = int(chimline.split("\t")[5])
                        chimfusendleft = int(chimline.split("\t")[6])
                        chimdepthleft = chimline.split("\t")[7]
                        chimfusstartright = chimline.split("\t")[9] 
                        chimfusendright = chimline.split("\t")[10]
                        chimdepthrigth = chimline.split("\t")[11]
                        if pseudogene == chimpseudogene:
                            if chimfuschrom == cfusionchrom:
                                if chimfusstartright != "NA" and int(chimfusstartleft) <= cfusionstart <= int(chimfusendright): # If we have both anchors the fusion should be inbetween the start left and the end right. 
                                    print >> out, chimline +"\t"+ cfusionchrom + "\t" + str(cfusionstart) + "\t" + str(cfusionend)  + "\t" + cfusiondepth
                                    clippedlist.append(repr(cline.rstrip())) # Append the entire row from both the chimread evidence and the clippedread evidence if we found them, if we did not detect them we can print them seperately as the final results 
                                    chimlist.append(repr(chimline.rstrip()))
                                    DetectedPseudogenesList.append(pseudogene)
                                elif chimfusstartright == "NA" and (int(chimfusstartleft) - chimreadpairdistance) <= cfusionstart <= (int(chimfusendleft) + chimreadpairdistance): # We dont have the right anchor and the clipp is within the left anchor or 100 plus or minus from the left anchor by default
                                    print >> out, chimline +"\t"+ cfusionchrom + "\t" + str(cfusionstart) + "\t" + str(cfusionend)  + "\t" + cfusiondepth
                                    clippedlist.append(repr(cline.rstrip()))
                                    chimlist.append(repr(chimline.rstrip()))
                                    DetectedPseudogenesList.append(pseudogene) 
        with open(PseudogeneCandidateChimbed, "r") as chimcoord:
            for chimline in chimcoord:
                #print repr(chimline.rstrip())
                #print chimlist
                rawstring = repr(chimline.rstrip())
                if not rawstring in chimlist: # In this case we have evidence from only a chimeric pair 
                    print >> out, chimline.rstrip() +"\tNA\tNA\tNA\tNA"
        with open(clippedwithpseudogeneoverlap, "r") as clippedCoord:
            for cline in clippedCoord:
                pseudogene = cline.split("\t")[0]
                pseudogenechrom = str(cline.split("\t")[1])
                pseudogenestart = str(cline.split("\t")[2])
                pseudogeneend = str(cline.split("\t")[3])
                cfusionchrom= str(cline.split("\t")[4])
                cfusionstart = str(cline.split("\t")[5])
                cfusionend = str(cline.split("\t")[6])
                cfusiondepth = str(cline.split("\t")[7])
                if not repr(cline.rstrip()) in clippedlist: # In this we have evidence from only clipped reads 
                    print >> out, pseudogene + "\t"+ pseudogenechrom + "\t" + pseudogenestart + "\t" + pseudogeneend + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t" + cfusionchrom + "\t" + cfusionstart + "\t" + cfusionend + "\t" + cfusiondepth.rstrip()
    return (clipandchimevidence, DetectedPseudogenesList)


def Pseuodogenecandidates(Sample,baminput,exoncoords, pseudogenecoords,chrominfo,genecoords,Cleaninglist, MovingList,Outputfolder,Psdepth):
    '''
    Extract reads that are split across exons, These are my pseudogene candidates  
    '''
    
    # Files

    Cigarbam = Sample + ".extractedN.bam"
    MovingList.append(Cigarbam)
    Cigarbai = Sample + ".extractedN.bam.bai"
    MovingList.append(Cigarbai)
    Cigarbed = Sample + ".extractedN.bed"
    Cleaninglist.append(Cigarbed)
    Pseudogenecandidates = Sample + "_pseudogenecandidates.tmp"
    Cleaninglist.append(Pseudogenecandidates)
    Pseudogenecandidatesbed = Sample + ".PseudogeneCandidates.bed"
    command_extractNcigar = "samtools view -h %s | awk '$1 ~ \"@\" || $6 ~ \"N\" {print $0}' | samtools view -Sb - | samtools sort -o %s" %(baminput, Cigarbam)
    os.system(command_extractNcigar)
    command_index = "samtools index %s" % Cigarbam
    os.system(command_index)
    logging.info('%s\tLooking for pseudogene Candidates', time.ctime().split(" ")[-2])
    logging.info('%s\tCalculating coverage over spliced alignment, filt 5 depth', time.ctime().split(" ")[-2])
    # Create a bed file by using genomecov for having a treshold in mapped reads! For obtaining pseudogene candidates
    commandgenomecovforsplittranscriptreads = "genomeCoverageBed -ibam %s -g %s -bg -split | awk '($4>%s){print $1\"\t\"$2\"\t\"$3}' > %s" %(Cigarbam, chrominfo,Psdepth,Cigarbed) # Atleast 5 in coverage 
    os.system(commandgenomecovforsplittranscriptreads) 
    logging.info('%s\tExtracting pseudogene candidates', time.ctime().split(" ")[-2])
    commandoverlappwithknowngenes = "intersectBed -wb -a %s -b %s -split | cut -f 7 | sort -n | uniq | cut -f 2- -d \"-\" | sort -n | uniq -c | sed -e 's/^[ \t]*//' | awk '$1 > 2 {print $2}' > %s" % (Cigarbed,exoncoords, Pseudogenecandidates) 
    os.system(commandoverlappwithknowngenes)
    if os.stat(Pseudogenecandidates).st_size==0: # If file is empty, no pseudogenes than you can end the script But first remove tmp files and then count the known processed pseudogenes 
        logging.info('%s\tCleaning...', time.ctime().split(" ")[-2])
        CountKnownPseudogenes(Sample,baminput,pseudogenecoords, MovingList)
        cleaning(Cleaninglist, MovingList,Outputfolder)
        sys.exit("No pseudogeneCandidates detected")
    else:
        Cleaninglist.append(Pseudogenecandidatesbed)
        logging.info('%s\tPseudogene candidates detected', time.ctime().split(" ")[-2])
        command = "grep -w -f %s %s > %s" %(Pseudogenecandidates, genecoords, Pseudogenecandidatesbed)
        os.system(command)
        return (Pseudogenecandidatesbed,Cigarbam)

def CountKnownPseudogenes(Sample,baminput,pseudogenecoords, MovingList):
    '''
    Here we annotate the known pseudogenes with bedtools 
    '''
    logging.info('%s\tCounting known pseudogenes... reads overlapping with the pseudogene coords in Ensembl',time.ctime().split(" ")[-2])
    AmountOfreadsinKnownProcessedPseudogenes = Sample + ".KnownProcessedPseudogenes.Out" 
    MovingList.append(AmountOfreadsinKnownProcessedPseudogenes)
     # Look at the intersect with known processed pseudogenes  
    bedintersectcommand = "intersectBed -abam %s -b %s -bed -wb" %(baminput, pseudogenecoords)
    capture_overlap_knownProcessedPseudogenes = subprocess.Popen(bedintersectcommand,shell = True, stdout=subprocess.PIPE)
    KnownProcessedPseudogenesList = [] # Create empty list for storing the pseudogene candidates 
    for line in capture_overlap_knownProcessedPseudogenes.stdout.readlines(): # Read each line of the standard out 
        gene = line.strip().split("\t")[-1] # Save last item in the output from intersect, the gene names for the overlapping reads 
        KnownProcessedPseudogenesList.append(gene) # Save to empty list 
    counts = Counter(KnownProcessedPseudogenesList) # Count the amount of times each gene appears, save it in a count object that works as a dictionary 
    
    min_treshold = 10 # Set treshold of how many reads the pseudogene most intersect with to be noted as there 
    filteredCounts = {x : counts[x] for x in counts if counts[x] >= min_treshold } # Filter the counts based on a treshold, amounts of reads that intersect with a known processed pseudogene!  
    # Get the positions of the pseudogene
    #filteredCountsSorted = sorted(filteredCounts, key=filteredCounts.get, reverse=True) # Sort the dictionary 
    DictWithPostAndAmountsOfReads = {}
    with open(pseudogenecoords, "r") as knownpseudogenes: # Open the pseudogene bed file
        for line in knownpseudogenes: 
            for key, value in filteredCounts.iteritems(): # Go through the dictionary 
                if line.strip().split("\t")[-1] == key:  # If gene name in known pseudogene file is the same as in the dictionary, append information about chr, start, end to the values
                    chrom = str(line.strip().split("\t")[0]) # append chromosome
                    Start = str(line.strip().split("\t")[1]) # append start
                    End = str(line.strip().split("\t")[2]) # append end
                    DictWithPostAndAmountsOfReads[key] = [chrom, Start, End, str(value)]
    with open(AmountOfreadsinKnownProcessedPseudogenes, "w") as out: 
        for key, values in DictWithPostAndAmountsOfReads.iteritems():
            out.write(key + "\t" + '\t'.join(values) + "\n")


def makeCircosGraph(Sample,baminput,Cigarbam,clipandchimevidence,exoncoords, chrominfo,DetectedPseudogenesList, Cleaninglist, MovingList): 
    '''
    Here we create the circos picture, Observe that i choose to create the circos pictures in the only case where the pseudogene is evident by both chimreads and pairs
    '''
    logging.info('%s\tCreating circos graph for samples with evidence from both Clipped reads and chimeric pairs...',time.ctime().split(" ")[-2])
    circosbin = "/apps/bio/apps/circos/0.67-7/bin/circos" # Exact path to circos     
    for i in DetectedPseudogenesList: #
        # Here we loop through our unique pseudogens for creating the intermediate files that are only required for the plotting of each pseudogene, exons and coverage for the ordinary bam and the split bam together with the link file, this will set up the backbone for the circos file 
        exactcoordbed = Sample + "_" + i + ".exoncoords.bed"
        Cleaninglist.append(exactcoordbed)
        covoverbam_circosformat = Sample + "_" + i + ".Cov_overBam.depth_circos"
        MovingList.append(covoverbam_circosformat)
        covoversplitbam_circosformat = Sample + "_" + i + ".Cov_overSplitBam.depth_circos"
        MovingList.append(covoversplitbam_circosformat) 
        karyotype = Sample + "_" + i + "_Karyotype"
        MovingList.append(karyotype)
        links = Sample + "_" + i + "_LinksbetweenExons.txt"
        MovingList.append(links)
        command = "grep -w %s$ %s > %s" %(i,exoncoords,exactcoordbed)
        os.system(command)

        # Creating one Karyotype for each Pseudogene structure and one link between each pseudogene
        with open(karyotype, "w") as karyotypeExons:
            with open(exactcoordbed, "r") as exoncoordread:
                for line in exoncoordread:
                    line = line.strip()
                    chromosome = str(line.split("\t")[0])
                    startexoncoord = int(line.split("\t")[1])
                    endexoncoord = int(line.split("\t")[2])
                    identifier = line.split("\t")[3]
                    print >> karyotypeExons, "chr - %s %s %s %s exon" %(identifier,identifier,startexoncoord-startexoncoord, abs(endexoncoord-startexoncoord))
        # Get the structure over the split Ordinary bam, immediately to circos format!
        command = "/home/xabras/Programs/bedtools2/bin/coverageBed -sorted -b %s -a %s -d -split -g %s| awk 'OFS=\" \" {print $4,$5-1,$5-1,$6}' > %s" %(baminput, exactcoordbed, chrominfo,covoverbam_circosformat) # OBS memory is consumed!!! That is why it breaks? 
        os.system(command)
        # Get the structure over the split N's bam, immediately to circos format!
        command = "/home/xabras/Programs/bedtools2/bin/coverageBed -sorted -b %s -a %s -d -split -g %s | awk 'OFS=\" \" {print $4,$5-1,$5-1,$6}' > %s" %(Cigarbam, exactcoordbed, chrominfo,covoversplitbam_circosformat)
        os.system(command)
        # Create the link file
        df_karyotype = pd.read_csv(karyotype, sep = " ", names=['chr', '-','identifier', 'identifier_2', 'start', 'end', 'color'])
        exonends = df_karyotype[df_karyotype['color'] == 'exon'][['identifier','end']][:-1] # get only exons from karyotype, then get only the identifier and end column and finally remove the last exons as you dont want a link from this one
        exonstarts = df_karyotype[df_karyotype['color'] == 'exon'][['identifier','start']][1:] # same as above but here we drop the first link because we dont want a start link from the first exon
        exonslinkdf = pd.concat([exonends.reset_index(drop=1).add_suffix('_1'), exonstarts.reset_index(drop=1).add_suffix('_2')], axis=1).dropna(axis=1,how="any")[['identifier_1','end_1', 'end_1','identifier_2','start_2','start_2']]
        exonslinkdf.to_csv(links, sep = " ", header =False, index = False)
        ## Create the config file
        configurationfile = Sample + "_" + i + ".conf"
        MovingList.append(configurationfile)
        outputpicture = Sample + "_" + i + ".png"
        outputpicturesvg = Sample + "_" + i + ".svg"
        Cleaninglist.append(outputpicture)
        Cleaninglist.append(outputpicturesvg)
        # Detect the max value from the regular bam coverage. This will set the treshold for the splits to follow
        depthlist = []
        with open(covoverbam_circosformat, "r") as d: # looping through the regular coverage from the bamfile over the exons. Save the largest value, this will be set as a max used both for the histogram and the heatmap
            for line in d: 
                line=line.strip()
                coveragebamdepth = line.split(" ")[3]
                depthlist.append(coveragebamdepth)
        maxcoverage = max(map(int,depthlist))
        # Here we print the configuration file
        with open(configurationfile, "w") as conf:
            print >> conf, """

# Plot the Circos figure 

# circos.conf

karyotype = %s

<<include colors_fonts_patterns.conf>>
<colors>
exon = 191,191,191

</colors>

<ideogram>
<spacing>
default = 0.005r
</spacing>
radius=0.6r
thickness=20p
fill=yes
show_label=yes # change here?
label_font=default
label_size=20p
label_radius=1.2r
label_parallel=no

</ideogram>

<image>
file*=%s
<<include etc/image.conf>>
</image>

# OBS you are switching the housekeeping conf here. This is due to all error messages
<<include /jumbo/WorkingDir/B17-006/PseudoScriptDb/housekeeping.conf>> 

# RGB/HSV color definition, color lists, location of fonts, fill patterns, included from circos distribution
<<include etc/colors_fonts_patterns.conf>>

show_ticks=yes
show_tick_labels=yes
<ticks>
radius=1r
color=black
thickness=5p
<tick>
show_labels=yes
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


            """ %(karyotype,outputpicture,links,covoversplitbam_circosformat,covoverbam_circosformat, maxcoverage) 
                
        # Plot the circos figures, here i use my path to circos  
        #command = "%s --conf %s " %(circosbin,configurationfile)
        command = "circos --conf %s" % configurationfile
        os.system(command)
        with open(clipandchimevidence, "r") as Sout: # Loop through the line of all the chimeric reads and ad the text on each and one of them. 
            next(Sout) # Skip the header
            for line in Sout: 
                line = line.strip() 
                SummaryPseudogene = line.split("\t")[0]
                chimreadsdepth = line.split("\t")[15]
                startFusion = line.split("\t")[5]
                chromFusion = line.split("\t")[4]
                CircosPictureWithAnno = Sample + "_" + line.replace("\t","_") + ".png"
                if SummaryPseudogene == i and line.split("\t")[5] != "NA" and line.split("\t")[13] != "NA":  # If the picture that we are looping over is the same as the pseudogene fusion, create the new png and move to output dict 
                    try: # If you have the left anchor add them together, otherwise just use the right anchor
                        chimpairsdepth = int(line.split("\t")[7]) + int(line.split("\t")[11])
                    except ValueError:
                        chimpairsdepth = int(line.split("\t")[7])
                    command = "convert -size 1500x1500 -gravity Center -pointsize 70 -annotate 0 \"%s -> %s:%s\\n\\n ChimReads:%s\\n\\n ChimPairs:%s\" %s %s" %(i,chromFusion,startFusion,chimreadsdepth,chimpairsdepth,outputpicture, CircosPictureWithAnno)
                    os.system(command)
                    MovingList.append(CircosPictureWithAnno)    
    

def GvizPlottingForOutput(Sample,baminput,Cigarbam,summaryOutAnnotatedBoth,exoncoords,Cleaninglist,MovingList):
    """
    This part plots the outputs using the gviz package in R instead of circos
    """
    logging.info('%s\tPlottingTheDetectedPseudogeneCandidates',time.ctime().split(" ")[-2])
    #bamfile=baminput
    #cigarbamfile=Cigarbam
    #outputReport=summaryOutAnnotatedBoth
    #exoncoordsdb = exoncoords
    with open(summaryOutAnnotatedBoth, "r") as finaloutputreport: 
        h=next(finaloutputreport)
        exoncoordsdict = {}
        for line in finaloutputreport:
            estarts = []
            eends = []
            line = line.strip()
            parentgene = line.split("\t")[0]
            parentchrom = line.split("\t")[1]
            parentgenestart = line.split("\t")[2]
            parentgeneend = line.split("\t")[3]
            if not line.split("\t")[4] == "NA":
                fuschrom  = line.split("\t")[4]
            else: # We only have evidence in the clip therefore you need to extract it from column 12
                fuschrom  = line.split("\t")[12]
            fusstartleft = line.split("\t")[5] 
            fusendleft = line.split("\t")[6]
            fusstartright = line.split("\t")[9]
            fusendright =  line.split("\t")[10]
            chimreadstart = line.split("\t")[13]
            chimreadend = line.split("\t")[14]
            Anno = line.split("\t")[16]+"_"+line.split("\t")[17] 
            command = "grep -w %s %s" %(parentgene,exoncoords)
            a=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            for excord in a.stdout:
                excord=excord.strip()
                exonchrom = excord.split("\t")[0]
                exonstart = int(excord.split("\t")[1])
                exonend = int(excord.split("\t")[2])
                estarts.append(exonstart)
                eends.append(exonend)
            if not parentgene in exoncoordsdict.keys(): 
                exoncoordsdict[parentgene]=(estarts,eends)
            else: 
                continue
            if not chimreadstart == "NA": 
                outputname = "%s_%swithin%s-%s" %(Sample,parentgene,fuschrom,chimreadstart)
                fusionstartplot = chimreadstart # When you have the clipp coord use them for the plotting
                fusionwidth = int(chimreadend)-int(chimreadstart)
                fusionrangestart = int(chimreadstart)-100
                fusionrangeend = int(chimreadstart)+100
            else: 
                outputname = "%s_%swithin%s-%s" %(Sample,parentgene,fuschrom,fusstartleft)
                if fusstartright == "NA": # You only have the left anchor, use those for the plotting
                    fusionstartplot = fusstartleft
                    fusionwidth = int(fusendleft) - int(fusstartleft)                
                    fusionrangestart = int(fusstartleft)-100 
                    fusionrangeend = int(fusendleft)+100
                else: # You have the entire anchor rang , use that for the plotting 
                    fusionstartplot = fusstartleft
                    fusionwidth = int(fusendright)-int(fusstartleft)
                    fusionrangestart = int(fusstartleft)-100
                    fusionrangeend = int(fusendright)+100
            outputRscript = outputname + ".R"
            outputPicturepdf = outputname + ".pdf" 
            outputPicturepng = outputname + ".png" # png for having in the ppsy reports
            with open(outputRscript, "w") as R:             
                startvector= "starts=c("+','.join(str(e) for e in exoncoordsdict[parentgene][0]) + ")"
                endvector= "ends=c("+','.join(str(e) for e in exoncoordsdict[parentgene][1]) + ")"          
                print >> R, """
#!/usr/bin/R  
suppressMessages(library(Gviz))
options(ucscChromosomeNames=FALSE)
axisTrack <- GenomeAxisTrack()              
%s 
%s
widths=ends-starts
bamfile <- '%s'
cigarbamfile <- '%s'
# AlignmentTracks 
alTrack <- DataTrack(range=bamfile, genome = "hg19", name = "Raw Coverage",  type="heatmap",chromosome="%s")
alTrack2 <- AlignmentsTrack(cigarbamfile,chromosome="%s",genome="hg19",name="Splice Coverage", isPaired=TRUE)
bamzoominsert<-DataTrack(range=bamfile, genome = "hg19", name = "Insert Coverage",  type="histogram", chromosome="%s")
# AnnotationTracks 
aTrack.groups <- AnnotationTrack(start=c(starts), width=c(widths),chromosome="%s", strand=rep("*",length(starts)),group=rep("%s",length(starts)), genome = "hg19", name = "Anno")
AnnotationTrackinsert<-AnnotationTrack(start=%s, width = %s, strand="*",shape="box",id="%s",genome="hg19", chromosome="%s",name ="Insert", group = "%s")

pdf("%s") 
grid.newpage()
pushViewport(viewport(layout=grid.layout(4, 6)))
pushViewport(viewport(layout.pos.col=c(1,2,3,4,5,6), layout.pos.row=c(1,2)))
plotTracks(list(axisTrack,alTrack,alTrack2,aTrack.groups), groupAnnotation="group", just.group=\"above\",from=%s,to=%s,type=c("heatmap","coverage"), sizes=c(1,2,3,0.5),fill.coverage=\"grey\",add=TRUE)
popViewport(1) 
pushViewport(viewport(layout.pos.col=c(2.5,5.5), layout.pos.row=4))
plotTracks(c(bamzoominsert,AnnotationTrackinsert,axisTrack), from=%s,to=%s, fill="darkred", sizes =c(0.5,0.1,0.5), lwd = 1, add=TRUE, panel.only=TRUE, legend=TRUE, groupAnnotation = "id", fill.histogram="darkgrey")    
dev.off()

png("%s") 
grid.newpage()
pushViewport(viewport(layout=grid.layout(4, 6)))
pushViewport(viewport(layout.pos.col=c(1,2,3,4,5,6), layout.pos.row=c(1,2)))
plotTracks(list(axisTrack,alTrack,alTrack2,aTrack.groups), groupAnnotation="group", just.group=\"above\",from=%s,to=%s,type=c("heatmap","coverage"), sizes=c(1,2,3,0.5),fill.coverage=\"grey\",add=TRUE)
popViewport(1) 
pushViewport(viewport(layout.pos.col=c(2.5,5.5), layout.pos.row=4))
plotTracks(c(bamzoominsert,AnnotationTrackinsert,axisTrack), from=%s,to=%s, fill="darkred", sizes =c(0.5,0.1,0.5), lwd = 1, add=TRUE, panel.only=TRUE, legend=TRUE, groupAnnotation = "id", fill.histogram="darkgrey")    
dev.off()


            """ %(startvector,endvector, baminput,Cigarbam ,parentchrom, parentchrom,fuschrom,parentchrom,parentgene,fusionstartplot,fusionwidth,Anno,fuschrom,Anno ,outputPicturepdf,parentgenestart,parentgeneend,fusionrangestart,fusionrangeend,outputPicturepng,parentgenestart,parentgeneend,fusionrangestart,fusionrangeend)

            command = "Rscript %s" %outputRscript # This part plots using the Gviz plotting
            os.system(command)
            # Cleaning up the outputs
            MovingList.append(outputPicturepdf)
            MovingList.append(outputPicturepng)
            MovingList.append(outputRscript)
            
def AnnotateFusionPointWithAnnovar(Sample,clipandchimevidence,anndb,annovarscript,Cleaninglist,MovingList): 
    '''
    Here we annotate the fusion point with annovar, using the absolute start position for the fusion range. 
    '''
    #filestoremove = []
    logging.info('%s\tAnnotating the fusion points using annovar',time.ctime().split(" ")[-2])
    summaryOutAnnotatedBoth = Sample + ".ChimPairs_ChimReads.Ppsy.txt"
    MovingList.append(summaryOutAnnotatedBoth)
    # Calculate the insert within a gene for the fusion 
    with open(summaryOutAnnotatedBoth, "w") as copyexactcoord:
        print >> copyexactcoord, "Pseudogene\tP_chrom\tP_start\tP_end\tFus_chrom_left\tFus_start_left\tFus_end_left\tnChimPairs_left\tFus_chrom_right\tFus_start_right\tFus_end_right\tnChimPairs_right\tChimRead_chrom\tChimRead_Start\tChimRead_End\tnChimReads\tInsertDescr\tInsertDescr"
        with open(clipandchimevidence, 'r') as exactcoord: # Loop through a copy of the input file.
            next(exactcoord)
            for line in exactcoord:
                stripped=line.strip()
                if not stripped.split("\t")[13] == "NA": # If we have the clipped fusion point that one sets the coords, otherwise uer use the start coord on the right anchor
                    fuschr = stripped.split("\t")[12]
                    fusstart = stripped.split("\t")[13]
                    fusend = stripped.split("\t")[14]
                else: 
                    fuschr = stripped.split("\t")[4]
                    fusstart = stripped.split("\t")[5]
                    fusend = stripped.split("\t")[6]
                outtmp=Sample+"_"+str(fuschr)+"_"+str(fusstart)+".tmp"
                with open(outtmp, "w") as out:
                    print >> out, "%s\t%s\t%s\t0\t0" %(fuschr, fusstart,fusend)
                commandRunAnnovar = "perl %s --geneanno --buildver hg19 %s %s -out %s.Annovar.txt" %(annovarscript, outtmp, anndb, outtmp)
                os.system(commandRunAnnovar)
                with open(outtmp+".Annovar.txt.variant_function", "r") as variantfunction:
                    for lineannovar in variantfunction:
                        strippedannovar = lineannovar.strip()
                        info = strippedannovar.split("\t")[0]
                        ExactAnno = strippedannovar.split("\t")[1]
                        print >> copyexactcoord, stripped + "\t" + info + "\t" + ExactAnno
                for filename in glob.glob(outtmp+"*"):
                    os.remove(filename)    
    return summaryOutAnnotatedBoth

def cleaning(Cleaninglist, MovingList,Outputfolder):
    # Itererating for files to be removed
    logging.info('%s\tCleaning',time.ctime().split(" ")[-2])
    if Cleaninglist:
        for item in list(set(Cleaninglist)):
            os.remove(item)
    # move all the files to this folder, then create the subfolders and move the files to that one.
    # If there is duplicates in the list, set to unique, it is ugly but some are added within a for loop  
    if MovingList:
        for item in list(set(MovingList)):
            os.rename(item, Outputfolder+"/"+item)
    # Now create the subdirectories 
    alignmentout = "%s/Alignments" %Outputfolder
    os.makedirs(alignmentout)
    # Alignments
    alignmentfiles = [f for f in os.listdir(Outputfolder) if f.endswith((".bam",".bai"))]
    if alignmentfiles:
        for f in alignmentfiles:
            try: 
                os.rename(Outputfolder+"/"+f, alignmentout+"/"+f)
            except OSError: 
                pass 
    # Gviz 
    GvizOut  = "%s/Plotting" %Outputfolder
    os.makedirs(GvizOut)
    Rscriptout = GvizOut + "/Scripts"
    os.makedirs(Rscriptout)
    covplots =  [f for f in os.listdir(Outputfolder) if f.endswith((".png", ".pdf"))]
    if covplots: 
        for f in covplots:
            try:
                os.rename(Outputfolder+"/"+f, GvizOut+"/"+f)
            except OSError:
                pass
    Rscripts =  [f for f in os.listdir(Outputfolder) if f.endswith(".R")]
    if Rscripts: 
        for f in Rscripts: 
            try: 
                os.rename(Outputfolder+"/"+f, Rscriptout+"/"+f)
            except OSError: 
                pass
    # Circos
    Circosout = "%s/Circos" %Outputfolder
    os.makedirs(Circosout)
    Circosout_conf = Circosout + "/Configs"
    os.makedirs(Circosout_conf)
    Circosout_pics = Circosout + "/CircosPictures" 
    os.makedirs(Circosout_pics)
    configfiles = [f for f in os.listdir(Outputfolder) if f.endswith(("_circos","_Karyotype", ".conf", "_LinksbetweenExons.txt"))]
    if configfiles: 
        for f in configfiles:
            try: 
                os.rename(Outputfolder+"/"+f, Circosout_conf+"/"+f)
            except OSError: 
                pass 
    circospics = [f for f in os.listdir(Outputfolder) if f.endswith(".png")]
    if circospics: # If you have something in the list
        for f in circospics:
            try: 
                os.rename(Outputfolder+"/"+f, Circosout_pics+"/"+f)
            except OSError: 
                pass 
    # OutPut Results  
    Logout = "%s/PpsyReports" %Outputfolder
    os.makedirs(Logout)
    logfiles = [f for f in os.listdir(Outputfolder) if f.endswith((".ChimPairs_ChimReads.Ppsy.txt", "_PsudogenewithFusUnfiltered.txt","BinnedPseudoFusExonOverlap.txt",".KnownProcessedPseudogenes.Out"))]
    if logfiles: # If you have something in the list
        for f in logfiles: 
            try:
                os.rename(Outputfolder+"/"+f, Logout+"/"+f)
            except OSError: 
                pass 

def main(baminput,Sample, Psdepth,insdistance,ChimPairDepthTresh,ChimPairBinningTresh,ChimReadDepthTresh, ChimReadBinningTresh,chimreadpairdistance):
    # Create the logger
    logging.basicConfig(level=logging.INFO)
    logging.info('%s\tStarting Ppsyfinder', time.ctime())
    (chrominfo,genecoords,pseudogenecoords,exoncoords,anndb,annovarscript,Cleaninglist,MovingList,Outputfolder)=database(baminput,Sample)
    (Pseudogenecandidatesbed,Cigarbam)=Pseuodogenecandidates(Sample,baminput,exoncoords,pseudogenecoords,chrominfo,genecoords,Cleaninglist, MovingList,Outputfolder,Psdepth)
    bamChimericPairs=extractingchimericReads(Sample, baminput, Cleaninglist, MovingList,insdistance)
    PseudogeneCandidateChimbed=ChimericReadsOverlapwithPseudogeneCandidates(Sample, bamChimericPairs,Pseudogenecandidatesbed, Cleaninglist, MovingList,ChimPairDepthTresh,ChimPairBinningTresh)
    bamClipped=extractingClipped(Sample,baminput, Cleaninglist,MovingList) 
    clippedbed=ClippedBamtoCoord(bamClipped,Cleaninglist, MovingList)
    clippedmapping=clippedreadbinning(clippedbed,Cleaninglist, MovingList,ChimReadDepthTresh,ChimReadBinningTresh)
    clippedwithpseudogeneoverlap=intersectClippedPseudogeneCandidates(Sample, clippedmapping, Pseudogenecandidatesbed, Cleaninglist, MovingList)
    (clipandchimevidence, DetectedPseudogenesList)=CombiningClippWithChimericReads(Sample, clippedwithpseudogeneoverlap, PseudogeneCandidateChimbed, Cleaninglist, MovingList,chimreadpairdistance)
    CountKnownPseudogenes(Sample,baminput,pseudogenecoords, MovingList)
    summaryOutAnnotatedBoth=AnnotateFusionPointWithAnnovar(Sample,clipandchimevidence,anndb,annovarscript,Cleaninglist,MovingList)
    GvizPlottingForOutput(Sample,baminput,Cigarbam,summaryOutAnnotatedBoth,exoncoords,Cleaninglist,MovingList)
    #makeCircosGraph(Sample,baminput,Cigarbam,clipandchimevidence,exoncoords,chrominfo,DetectedPseudogenesList, Cleaninglist, MovingList)
    cleaning(Cleaninglist, MovingList,Outputfolder)
    logging.info('%s\tPpsy Finished, results in output folder',time.ctime().split(" ")[-2])


if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.baminput, arguments.Sample, arguments.Psdepth, arguments.insdistance,arguments.ChimPairDepthTresh,arguments.ChimPairBinningTresh,arguments.ChimReadDepthTresh,arguments.ChimReadBinningTresh, arguments.chimreadpairdistance)
