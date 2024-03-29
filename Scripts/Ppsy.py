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
import psutil
import re 

from argparse import ArgumentParser, RawTextHelpFormatter


# Debugging, if we don't have any candidates i think we should still create the empty file about the pseudogene candidates? 

def parseArgs():
    parser = argparse.ArgumentParser(description='Detects processed pseudogenes by looking at DNA data that splits across the splice junctions',formatter_class=RawTextHelpFormatter)    
    choices_helper = { "Fastq": "Input is Fastq, The pipeline will run STAR for you, method specific parameters are [-R1 Path to Read1, -R2, Path to Read2, -STARindex Path to human reference genome]",
                       "Bam": "Input is Bam which means that you have already aligned your fastq files towards the reference genome in the required way (Read more about how to align using the method description), method specific parameter is [-I, Path to the bam input]"}
    parser.add_argument('--Method',dest="method", choices=choices_helper,help='\n'.join("{}: {}".format(key, value) for key, value in choices_helper.items()), required=True)
    parser.add_argument('-S', dest='Sample', help='Sample name, descides the prefix of your outputs together with the output folder (required)', required=True)
    parser.add_argument('--pseudoCandidateDepth', dest='Psdepth', help='The minimum depth that supports the splice junctions, these are the resulting processed pseudogene candidates (default 5)', default=5 ,type=int) 
    parser.add_argument('--InsertDistance', dest='insdistance', help='What is the distance from the parent gene where we can have an pseudogene, low distance might increase the amount of detected pseudogenes but will also increase the amount of false positives. A low distance and you might hit inserted pseudogenes in the parent gene itself which is not very likely (default 200 000)', default=200000,type=int)
    parser.add_argument('--ChimericPairDepthTreshold', dest='ChimPairDepthTresh', help='The minimum amount of reads to suppport the chimeric pairs in the left anchor, (min 5), (default 10)', default=10,type=int)
    parser.add_argument('--ChimericPairBinningTreshold', dest='ChimPairBinningTresh', help='When the chimeric pairs are binned into the anchors the binning distance defines the distance for the read to belong to same bin, (default 500)', default=500,type=int)
    parser.add_argument('--ChimericReadDepthTreshold', dest='ChimReadDepthTresh', help='The minimum amount of reads to support the chimeric reads in the fusion site, (min 5) (default 10)', default=10,type=int)     
    parser.add_argument('--ChimericReadBinningTreshold', dest='ChimReadBinningTresh', help='When the chimeric reads are binned into the anchors the binning distance defines the distance for the read to belong to same bin', default=10,type=int)
    parser.add_argument('--MergeChimReadWithChimpairTresh', dest='chimreadpairdistance', help='When we are combining the results from the chimeric pairs and the chimeric reads we combine them if the chimeric read are withing the chimeric pair anchors or the chimeric read are in a user defined distance from the left anchor (default 100)', default=100,type=int)
    opts, rem_args = parser.parse_known_args()
    if opts.method == "Bam":
        parser.add_argument('-I', dest='baminput',help='Bamfile containing alignment performed by a splice aware aligner, need index in the same folder (required)',required=True)
    else: 
        parser.add_argument('-R1', dest='Fastq1',help='R1 from a Fastq file (required)',required=True) 
        parser.add_argument('-R2', dest='Fastq2',help='R2 from a Fastq file (required)',required=True)
        parser.add_argument('-STARindex', dest='STARindex',help='Path to the star index of the reference genome',required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def CheckingPaths(): 
    """
    This is the part which checks so you have all the correct paths to the tools that you are using 
    """
    # Testing Bedtools 
    command = "which bedtools"
    a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    stdout, stderr = a.communicate()
    if not stdout:
        print("Error, there is no bedtools in your path, Please install! Exiting!")
        sys.exit()
    # Testing Samtools, Impportant that you are not using two old version of samtools 
    command = "which samtools"
    a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    stdout, stderr = a.communicate()
    if not stdout:
        print("Error, there is no samtools in your path, Please install! Exiting!")
        sys.exit()
    # Testing R 
    #command = "which R"
    #a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    #stdout, stderr = a.communicate()
    #if not stdout:
     #   print "Error, there is no R in your path, Please install! Exiting!"
     #   sys.exit()
    # Test required R packages 
    #command = "Rscript -e 'library(\"Gviz\")'"
    #a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    #stdout, stderr = a.communicate()
    #if not stdout:
    #    print "Error, there GVIZ is not installed in your R version, Please install! Exiting!"
    #    sys.exit()
    # Testing convert 
    command = "which convert"
    a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    stdout, stderr = a.communicate()
    if not stdout:
        print("Error, there is no convert in your path, Please install! Exiting!")
        sys.exit()
    # Testing awk 
    command = "which awk"
    a = subprocess.Popen(command, stdout=subprocess.PIPE,stderr = subprocess.PIPE, shell=True)
    stdout, stderr = a.communicate()
    if not stdout:
        print("Error, there is no awk in your path, Please install! Exiting!")
        sys.exit()

def database(): 
    """
    This part reads the databases that is used for input 
    """
    # Check database after installation 
    logging.info('%s\tCollecting the paths to the HG19 databases', time.ctime().split(" ")[-2])
    genecoords = os.path.abspath(__file__).split("bin/Ppsy.py")[0] + "HG19_databases/Gene_coord_hg19_refgene.bed"
    exoncoords = os.path.abspath(__file__).split("bin/Ppsy.py")[0] + "HG19_databases/Exon_coord_hg19_refgene.bed"
    #genecoords = "/home/xabras/.conda/envs/Ppsy/HG19_databases/Gene_coord_hg19_refgene.bed"
    #exoncoords = "/home/xabras/.conda/envs/Ppsy/HG19_databases/Exon_coord_hg19_refgene.bed"
    pseudogenecoords = os.path.abspath(__file__).split("bin/Ppsy.py")[0] + "HG19_databases/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed"
    #pseudogenecoords="/home/xabras/.conda/envs/Ppsy/HG19_databases/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed"
    return(genecoords,pseudogenecoords,exoncoords)

def CreatingOutputDir(Sample):
    """
    This Part creates the output folder together with the cleaning lists 
    """
    logging.info('%s\tSetting Up Outputdir', time.ctime().split(" ")[-2])
    Outputfolder=Sample + "_PPsyOut"
    command = "mkdir -p %s" % Outputfolder
    os.system(command)
    # List for cleaning and moving files
    Cleaninglist=[]
    MovingList = []
    return(Cleaninglist, MovingList,Outputfolder)

def AlignmentWithSTAR(Fastq1,Fastq2,STARindex,Sample,Outputfolder): 
    """
    This function is only run if you dont have the alignment file so this part runs STAR, the Issue here is that we cannot run this in multicore while the other is single core, this is why it is preferable to have the Alignmentfile with star already done  
    """
    # You need to have STAR in your path for this 
    logging.info('%s\tAlignment using STAR... This step might be very slow', time.ctime().split(" ")[-2])
    baminput = Outputfolder +"/"+ Sample + "Aligned.sortedByCoord.out.bam"
    if Fastq1.endswith(".gz"):
        command = "STAR --runThreadN 1 --genomeDir %s --chimOutType WithinBAM --outSAMunmapped Within --outFilterMultimapNmax 20 --chimSegmentMin 20 --readFilesCommand zcat --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s/%s" %(STARindex,Fastq1,Fastq2,Outputfolder,Sample)
        os.system(command)
    else:
        command = "STAR --runThreadN 1 --genomeDir %s --chimOutType WithinBAM --outSAMunmapped Within --outFilterMultimapNmax 20 --chimSegmentMin 20 --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s/%s" %(STARindex,Fastq1,Fastq2,Outputfolder,Sample)
        os.system(command)
    command = "samtools index %s" %baminput
    os.system(command)
    return(baminput)

def CreateTheChromosomeInfoFile(baminput,Sample,Cleaninglist):
    """
    This extract the chromsome info file from the input 
    """
    logging.info('%s\tCreating Chromsome Info File', time.ctime().split(" ")[-2])
    chrominfo = Sample + ".ChromInfo_SortingOrder.txt"    
    command = "samtools view %s -H" % baminput 
    a = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    with open(chrominfo, "w") as chromout: 
        for line in a.stdout:
            line = line.strip().decode("utf-8")
            if line.split("\t")[0] == "@SQ":
                chrom  = line.split("\t")[1].split(":")[-1]
                size = line.split("\t")[2].split(":")[-1]
                print(chrom +"\t"+ str(size), file=chromout)
    Cleaninglist.append(chrominfo)
    return chrominfo
    

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
                            print(chrom + "\t" + str(startInGene) + "\t" + str(endInGene)+"\t"+matechrom + "\t" + str(startMate) + "\t" + str(endMate), file=tmpout)
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
                    sorted[["diff1", "diff2"]] = sorted.groupby(3)[[1,4]].diff()
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
                                                    print(pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\t"+str(row.loc[3])+ "\t" + str(FusionAnchorrightStart)+"\t" + str(FusionAnchorrightEnd) + "\t" + str(len(PseudoanchorRightstarts)), file=out)
                                                    
                                                else: 
                                                    # We did not detect the second anchor
                                                    print(pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\tNA\tNA\tNA\tNA", file=out)
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
                                        leftstring = str(row.loc[3]) + "\t" + str(FusionAnchorleftStart) + "\t" + str(FusionAnchorleftEnd) + "\t"+str(len(FusionAnchorleftstarts))
                                        if len(PseudoanchorRightstarts) > ChimPairDepthTresh:
                                            # This part tells us we have the right anchor as well 
                                            PseudogeneAnchorrightStart=PseudoanchorRightstarts[0]
                                            PseudogeneAnchorrightEnd=PseudoanchorRightends[-1]
                                            FusionAnchorrightStart=FusionAnchorRightstarts[0]
                                            FusionAnchorrightEnd=FusionAnchorRightends[-1]
                                            print (pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\t"+str(row.loc[3])+ "\t" + str(FusionAnchorrightStart)+"\t" + str(FusionAnchorrightEnd) + "\t" + str(len(PseudoanchorRightstarts)), file=out)
                                        
                                        else: 
                                            # We did not detect the second anchor
                                            print(pseudogenecandidate+"\t"+ chrom + "\t" + str(start) + "\t" + str(end)+ "\t"+leftstring + "\tNA\tNA\tNA\tNA", file=out)
                                    
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
                    Cigar = line.split("\t")[5]
                    #print(list(Cigar))
                    if not "D" in Cigar or not "I" in Cigar: # We are ignoring the reads containing a D or an I 
                        #print(Cigar)
                        char= re.sub('\d+', ',', Cigar)
                        char=char.split(",")
                        char=list(filter(None, char))
                        Cigar=re.sub('\D', ',', Cigar)
                        Cigar=Cigar.split(",")
                        Cigar = list(filter(None, Cigar))
                        CigarList = [int(i) for i in Cigar] 
                        startchimCigar= line.split("\t")[16].split(":")[2].split(",")[3]
                        charchim=re.sub('\d+', ',', startchimCigar)
                        charchim=charchim.split(",")
                        charchim=list(filter(None, charchim))
                        cigarchim=re.sub('\D', ',', startchimCigar)
                        cigarchim=cigarchim.split(",")
                        cigarchim=list(filter(None, cigarchim))
                        cigarchim=[int(i) for i in cigarchim]
                        #print(charchim)
                        #print(cigarchim)
                        if len(cigarchim)==3: # If we have longer then 2 remove the softclipped as you wont count it anyway 
                            minval=cigarchim.index(min(cigarchim))
                            charchim.pop(minval)
                            cigarchim.pop(minval)
                        if len(CigarList)==3:
                            minval=CigarList.index(min(CigarList))
                            CigarList.pop(minval)
                            char.pop(minval)
                        if len(cigarchim) == 2 and len(CigarList)==2: # If both of them are 2 it is very straight forward! 
                            startchim = line.split("\t")[16].split(":")[2].split(",")[1]
                            if "M" == charchim[0]: 
                                startchim=str(int(startchim)+cigarchim[0])
                            else:
                                startchim=str(int(startchim))
                            if "M" == char[0]:
                                start=str(int(line.split("\t")[3])+int(CigarList[0]))
                                print(chrom+":"+start+"\t"+chromchim+":"+startchim, file=out)
                                #print(chrom+":"+start+"\t"+chromchim+":"+startchim)
                                listofreads.append(read)
                            else:
                                start=str(int(line.split("\t")[3]))
                                print(chrom+":"+start+"\t"+chromchim+":"+startchim, file=out)
                                #print(chrom+":"+start+"\t"+chromchim+":"+startchim)
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
                    for key, value in d.items():
                        firstitem = 0 
                        amount = 1 
                        if len(value) > ChimReadDepthTresh: # If you have atleast 10 in depth by default 
                            sortedlist = sorted(map(int,value))
                            for i in range(1,len(sortedlist)): # Remeber the indexes
                                if sortedlist[i] - sortedlist[i-1] < ChimReadBinningTresh: # If the previus item is atleast 20 away from the new one by default
                                    amount += 1 
                                    if i+1 == len(value) and amount > ChimReadDepthTresh: # If you are reaching the end of the list and you have atleast 5 
                                        print(line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i]), file=clippedout)
                                        #print(line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i]))
                                else: 
                                    if amount > ChimReadDepthTresh: # not in the end of the list but you have 10 by default
                                        print(line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i-1]),file=clippedout)
                                        #print(line + "\t" + str(amount) + "\t" + key + "\t" + str(sortedlist[firstitem]) + "\t" + str(sortedlist[i-1]))
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
                        #print(cline)
                        clipfirstchrom = cline.split("\t")[0] 
                        clipsecondchrom = cline.split("\t")[4] 
                        clipfirststart = int(cline.split("\t")[1] )
                        clipsecondstart = int(cline.split("\t")[5] )
                        clipdepth = cline.split("\t")[3] 
                        clipfirstend = int(cline.split("\t")[2] )
                        clipsecondend = int(cline.split("\t")[6])
                        if pseudochrom == clipfirstchrom: # If first coordinates in the clip overlaps the pseudogene candidate 
                            if pseudostart-10 < clipfirststart < pseudoend+10 or pseudostart-10 < clipfirstend < pseudoend+10: 
                                print(pseudogene+"\t"+str(pseudochrom)+"\t"+str(pseudostart)+"\t"+str(pseudoend)+"\t"+str(clipsecondchrom)+"\t"+str(clipsecondstart)+"\t"+str(clipsecondend)+"\t"+str(clipdepth),file=out)
                        elif pseudochrom == clipsecondchrom:       
                            if pseudostart-10 < clipsecondstart < pseudoend+10 or pseudostart-10 < clipsecondend < pseudoend+10:
                                print(pseudogene+"\t"+str(pseudochrom)+"\t"+str(pseudostart)+"\t"+str(pseudoend)+"\t"+str(clipfirstchrom)+"\t"+str(clipfirststart)+"\t"+str(clipfirstend)+"\t"+str(clipdepth), file=out)
    if not os.stat(clippedwithpseudogeneoverlap).st_size==0: # The clipped reads with the pseudogenes might be empty and this will break the loop. Therefore we are making sure that the file is not empty before the analysis
        pandasdf=pd.read_csv(clippedwithpseudogeneoverlap, sep = "\t", header=None) # If we have the same fusion in different positions in the pseudogene we merge them 
        pandasdf2= pandasdf.groupby([0,1,2,3,4,5,6], as_index=False)[7].sum()
        #print(pandasdf2)
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
        print("Pseudogene\tP_chrom\tP_start\tP_end\tIS_chrom_left\tIS_start_left\tIS_end_left\tnChimPairs_left\tIS_chrom_right\tIS_start_right\tIS_end_right\tnChimPairs_right\tChimRead_chrom\tChimRead_Start\tChimRead_End\tnChimReads", file=out)
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
                                    print(chimline +"\t"+ cfusionchrom + "\t" + str(cfusionstart) + "\t" + str(cfusionend)  + "\t" + cfusiondepth, file=out)
                                    clippedlist.append(repr(cline.rstrip())) # Append the entire row from both the chimread evidence and the clippedread evidence if we found them, if we did not detect them we can print them seperately as the final results 
                                    chimlist.append(repr(chimline.rstrip()))
                                    DetectedPseudogenesList.append(pseudogene)
                                elif chimfusstartright == "NA" and (int(chimfusstartleft) - chimreadpairdistance) <= cfusionstart <= (int(chimfusendleft) + chimreadpairdistance): # We dont have the right anchor and the clipp is within the left anchor or 100 plus or minus from the left anchor by default
                                    print(chimline +"\t"+ cfusionchrom + "\t" + str(cfusionstart) + "\t" + str(cfusionend)  + "\t" + cfusiondepth, file=out)
                                    clippedlist.append(repr(cline.rstrip()))
                                    chimlist.append(repr(chimline.rstrip()))
                                    DetectedPseudogenesList.append(pseudogene) 
        with open(PseudogeneCandidateChimbed, "r") as chimcoord:
            for chimline in chimcoord:
                #print repr(chimline.rstrip())
                #print chimlist
                rawstring = repr(chimline.rstrip())
                if not rawstring in chimlist: # In this case we have evidence from only a chimeric pair 
                    print(chimline.rstrip() +"\tNA\tNA\tNA\tNA", file=out)
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
                    print(pseudogene + "\t"+ pseudogenechrom + "\t" + pseudogenestart + "\t" + pseudogeneend + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t" + cfusionchrom + "\t" + cfusionstart + "\t" + cfusionend + "\t" + cfusiondepth.rstrip(), file=out)
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
        # I will still add the output for the pseudogene report but i will keep it empty except for the header, then it will be included in the make ppsy report 
        summaryOutAnnotatedBoth = Sample + ".ChimPairs_ChimReads.Ppsy.txt" # This is what the output is called
        MovingList.append(summaryOutAnnotatedBoth)
        with open(summaryOutAnnotatedBoth, "w") as copyexactcoord:
            print("Pseudogene\tP_chrom\tP_start\tP_end\tIS_chrom_left\tIS_start_left\tIS_end_left\tnChimPairs_left\tIS_chrom_right\tIS_start_right\tIS_end_right\tnChimPairs_right\tChimRead_chrom\tChimRead_Start\tChimRead_End\tnChimReads\tIS_Region\tIS_Gene", file=copyexactcoord)
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
    AmountOfreadsinKnownProcessedPseudogenes = Sample + ".KnownProcessedPseudogenes.txt" 
    MovingList.append(AmountOfreadsinKnownProcessedPseudogenes)
    # Look at the intersect with known processed pseudogenes  
    #print pseudogenecoords
    bedintersectcommand = "intersectBed -abam %s -b %s -bed -wb" %(baminput, pseudogenecoords)
    capture_overlap_knownProcessedPseudogenes = subprocess.Popen(bedintersectcommand,shell = True, stdout=subprocess.PIPE)
    KnownProcessedPseudogenesList = [] # Create empty list for storing the pseudogene candidates 
    for line in capture_overlap_knownProcessedPseudogenes.stdout.readlines(): # Read each line of the standard out 
        gene = line.strip().decode("utf-8").split("\t")[-1] # Save last item in the output from intersect, the gene names for the overlapping reads 
        KnownProcessedPseudogenesList.append(gene) # Save to empty list 
    counts = Counter(KnownProcessedPseudogenesList) # Count the amount of times each gene appears, save it in a count object that works as a dictionary 
    
    min_treshold = 10 # Set treshold of how many reads the pseudogene most intersect with to be noted as there 
    filteredCounts = {x : counts[x] for x in counts if counts[x] >= min_treshold } # Filter the counts based on a treshold, amounts of reads that intersect with a known processed pseudogene!  
    # Get the positions of the pseudogene
    #filteredCountsSorted = sorted(filteredCounts, key=filteredCounts.get, reverse=True) # Sort the dictionary 
    DictWithPostAndAmountsOfReads = {}
    with open(pseudogenecoords, "r") as knownpseudogenes: # Open the pseudogene bed file
        for line in knownpseudogenes: 
            for key, value in filteredCounts.items(): # Go through the dictionary 
                if line.strip().split("\t")[-1] == key:  # If gene name in known pseudogene file is the same as in the dictionary, append information about chr, start, end to the values
                    chrom = str(line.strip().split("\t")[0]) # append chromosome
                    Start = str(line.strip().split("\t")[1]) # append start
                    End = str(line.strip().split("\t")[2]) # append end
                    DictWithPostAndAmountsOfReads[key] = [chrom, Start, End, str(value)]
    with open(AmountOfreadsinKnownProcessedPseudogenes, "w") as out: 
        for key, values in DictWithPostAndAmountsOfReads.items():
            out.write(key + "\t" + '\t'.join(values) + "\n")

def GvizPlottingForOutput(Sample,baminput,Cigarbam,summaryOutAnnotatedBoth,exoncoords,Cleaninglist,MovingList):
    """
    This part plots the outputs using the gviz package in R instead of circos
    """
    logging.info('%s\tPlotting The Detected PseudogeneCandidates using GVIZ',time.ctime().split(" ")[-2])
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
                excord=excord.strip().decode("utf-8")
                exonchrom = excord.split("\t")[0]
                exonstart = int(excord.split("\t")[1])
                exonend = int(excord.split("\t")[2])
                estarts.append(exonstart)
                eends.append(exonend)
            if not parentgene in exoncoordsdict.keys(): 
                exoncoordsdict[parentgene]=(estarts,eends)
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
            outputPicturepng = outputname + ".png" # png for having in the ppsy reports, cannot create them at the nodes so you need to do it with a seperate command
            with open(outputRscript, "w") as R:             
                startvector= "starts=c("+','.join(str(e) for e in exoncoordsdict[parentgene][0]) + ")"
                endvector= "ends=c("+','.join(str(e) for e in exoncoordsdict[parentgene][1]) + ")"
                print("""
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

TextTrack1 <- CustomTrack(plottingFunction = function(GdObject, prepare) {if(!prepare) grid.text("Coverage over the insert site") ; return(invisible(GdObject))}, name=NULL,displayPars=list(cex=1.5))

displayPars(aTrack.groups)<-list(fill="black", col=NULL)
displayPars(AnnotationTrackinsert)<-list(col=NULL)

pdf("%s") 
grid.newpage()
pushViewport(viewport(layout=grid.layout(4, 6)))
pushViewport(viewport(layout.pos.col=c(1,2,3,4,5,6), layout.pos.row=c(1,2)))
plotTracks(list(axisTrack,alTrack,alTrack2,aTrack.groups), groupAnnotation="group", just.group=\"above\",from=%s,to=%s,type=c("heatmap","coverage"), sizes=c(1,2,3,0.5),fill.coverage=\"grey\",add=TRUE, main="Coverage over the Gene",background.title="white",col.title="black",col.axis="black", cex.main=1.5)
popViewport(1) 
pushViewport(viewport(layout.pos.col=c(2.5,5.5), layout.pos.row=c(3,4)))
plotTracks(c(TextTrack1,bamzoominsert,AnnotationTrackinsert,axisTrack), from=%s,to=%s, fill="darkred", sizes =c(0.5,0.5,0.1,0.5), lwd = 1, add=TRUE, background.title="white",col.title="black",col.axis="black", groupAnnotation = "id", fill.histogram="darkgrey",col.histogram="darkgrey")
                      shhh<-dev.off()""" %(str(startvector),str(endvector), str(baminput),str(Cigarbam) ,str(parentchrom), str(parentchrom),str(fuschrom),str(parentchrom),str(parentgene),str(fusionstartplot),str(fusionwidth),str(Anno),str(fuschrom),str(Anno) ,str(outputPicturepdf),str(parentgenestart),str(parentgeneend),str(fusionrangestart),str(fusionrangeend)), file=R)
            command = "Rscript %s" %outputRscript # This part plots using the Gviz plotting
            os.system(command)
            # Cleaning up the outputs
            MovingList.append(outputPicturepdf)
            # Creating png does not work on the nodes, use the command convert instead for those ... 
            command = "convert %s %s" %(outputPicturepdf, outputPicturepng)
            os.system(command)
            MovingList.append(outputPicturepng)
            MovingList.append(outputRscript)
            
def AnnotateFusionPoint(Sample,clipandchimevidence,exoncoords,genecoords,Cleaninglist,MovingList): 
    '''
    Here we annotate the fusion point with annovar, using the absolute start position for the fusion range. 
    '''
    #filestoremove = []
    logging.info('%s\tAnnotating the fusion points',time.ctime().split(" ")[-2])
    summaryOutAnnotatedBoth = Sample + ".ChimPairs_ChimReads.Ppsy.txt"
    MovingList.append(summaryOutAnnotatedBoth)
    # Calculate the insert within a gene for the fusion 
    with open(summaryOutAnnotatedBoth, "w") as copyexactcoord:
        print("Pseudogene\tP_chrom\tP_start\tP_end\tIS_chrom_left\tIS_start_left\tIS_end_left\tnChimPairs_left\tIS_chrom_right\tIS_start_right\tIS_end_right\tnChimPairs_right\tChimRead_chrom\tChimRead_Start\tChimRead_End\tnChimReads\tIS_Region\tIS_Genes", file=copyexactcoord)
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
                with open(genecoords, "r") as genecoordforinsertanno: # Here we use the gene coord first, if on gene coord is detected grep the exons with the same gene coord from the alignment 
                    for l in genecoordforinsertanno:
                        l = l.strip()
                        genechrom = l.split("\t")[0]
                        genestart = l.split("\t")[1]
                        geneend = l.split("\t")[2]
                        genename = l.split("\t")[3]
                        if fuschr == genechrom and int(fusstart) >= int(genestart) and int(fusend) <= int(geneend):
                            geneofinsert = genename
                            command = "grep %s %s" % (genename,exoncoords)
                            pipe = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE)
                            for exoncoordline in pipe.stdout: # If we have an insert in a gene we grep the exon coordinates 
                                exoncoordline = exoncoordline.strip().decode("utf-8")
                                exonchrom = exoncoordline.split("\t")[0] # I need to ad the exon chrom as an overlapp importance as well as i cannot grep exactly as the Gene coord is not exactly the same as the exon coords which have an id onto it. 
                                exonstart = exoncoordline.split("\t")[1]
                                exonend = exoncoordline.split("\t")[2]
                                if fuschr == exonchrom and int(fusstart) >= int(exonstart) and int(fusend) <= int(exonend):
                                    geneinsertanno = "exonic"
                                    break 
                                else: 
                                    geneinsertanno = "intronic"
                            break 
                        else: # No hit within the genecoord therefore it is intergenic
                            geneofinsert = "NA"
                            geneinsertanno = "intergenic"  
                    print(stripped + "\t" + geneinsertanno + "\t" + geneofinsert,  file=copyexactcoord)
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
    if os.path.exists(alignmentout):
        shutil.rmtree(alignmentout)
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
    if os.path.exists(GvizOut):
        shutil.rmtree(GvizOut)
    os.makedirs(GvizOut)
    Rscriptout = GvizOut + "/Scripts"
    if os.path.exists(Rscriptout):
        shutil.rmtree(Rscriptout)
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
    # OutPut Results  
    Logout = "%s/PpsyReports" %Outputfolder
    if os.path.exists(Logout):
        shutil.rmtree(Logout)
    os.makedirs(Logout)
    logfiles = [f for f in os.listdir(Outputfolder) if f.endswith((".ChimPairs_ChimReads.Ppsy.txt", "_PsudogenewithFusUnfiltered.txt","BinnedPseudoFusExonOverlap.txt",".KnownProcessedPseudogenes.txt"))]
    if logfiles: # If you have something in the list
        for f in logfiles: 
            try:
                os.rename(Outputfolder+"/"+f, Logout+"/"+f)
            except OSError: 
                pass 

def mainBam(baminput,Sample,Psdepth,insdistance,ChimPairDepthTresh,ChimPairBinningTresh,ChimReadDepthTresh, ChimReadBinningTresh,chimreadpairdistance):
    # Create the logger
    process = psutil.Process(os.getpid())
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('%s\tStarting Ppsyfinder', time.ctime())
    CommandLine  = "%s --Method Bam --pseudoCandidateDepth %s --InsertDistance %s --ChimericPairDepthTreshold %s --ChimericPairBinningTreshold %s --ChimericReadDepthTreshold %s --ChimericReadBinningTreshold %s --MergeChimReadWithChimpairTresh %s -S %s -I %s" %(sys.argv[0], arguments.Psdepth,arguments.insdistance, arguments.ChimPairDepthTresh, arguments.ChimPairBinningTresh, arguments.ChimReadDepthTresh, arguments.ChimReadBinningTresh, arguments.chimreadpairdistance, arguments.Sample, arguments.baminput)
    logging.info('%s\tCommandline\t%s', time.ctime().split(" ")[-2],CommandLine)
    CheckingPaths()
    (genecoords,pseudogenecoords,exoncoords)=database()
    (Cleaninglist, MovingList,Outputfolder)=CreatingOutputDir(Sample)
    chrominfo=CreateTheChromosomeInfoFile(baminput,Sample,Cleaninglist)
    (Pseudogenecandidatesbed,Cigarbam)=Pseuodogenecandidates(Sample,baminput,exoncoords,pseudogenecoords,chrominfo,genecoords,Cleaninglist, MovingList,Outputfolder,Psdepth)
    bamChimericPairs=extractingchimericReads(Sample, baminput, Cleaninglist, MovingList,insdistance)
    PseudogeneCandidateChimbed=ChimericReadsOverlapwithPseudogeneCandidates(Sample, bamChimericPairs,Pseudogenecandidatesbed, Cleaninglist, MovingList,ChimPairDepthTresh,ChimPairBinningTresh)
    bamClipped=extractingClipped(Sample,baminput, Cleaninglist,MovingList) 
    clippedbed=ClippedBamtoCoord(bamClipped,Cleaninglist, MovingList)
    clippedmapping=clippedreadbinning(clippedbed,Cleaninglist, MovingList,ChimReadDepthTresh,ChimReadBinningTresh)
    clippedwithpseudogeneoverlap=intersectClippedPseudogeneCandidates(Sample, clippedmapping, Pseudogenecandidatesbed, Cleaninglist, MovingList)
    (clipandchimevidence, DetectedPseudogenesList)=CombiningClippWithChimericReads(Sample, clippedwithpseudogeneoverlap, PseudogeneCandidateChimbed, Cleaninglist, MovingList,chimreadpairdistance)
    CountKnownPseudogenes(Sample,baminput,pseudogenecoords, MovingList)
    summaryOutAnnotatedBoth=AnnotateFusionPoint(Sample,clipandchimevidence,exoncoords,genecoords,Cleaninglist,MovingList)
    GvizPlottingForOutput(Sample,baminput,Cigarbam,summaryOutAnnotatedBoth,exoncoords,Cleaninglist,MovingList)
    #makeCircosGraph(Sample,baminput,Cigarbam,clipandchimevidence,exoncoords,chrominfo,DetectedPseudogenesList, Cleaninglist, MovingList)
    cleaning(Cleaninglist, MovingList,Outputfolder)
    logging.info('%s\tPpsy Finished, results in output folder',time.ctime().split(" ")[-2])
    logging.info('Total memory Consumed in Bytes:\t%s',(process.memory_info()[0]))
    print(process.memory_info())
    logging.info('RunTime in Sec:\t%s', int(time.time()-start))


def mainFastq(Fastq1,Fastq2,STARindex,Sample,Psdepth,insdistance,ChimPairDepthTresh,ChimPairBinningTresh,ChimReadDepthTresh, ChimReadBinningTresh,chimreadpairdistance):
    """
    This is the Main function to run when you have a Fastq as input and need to run STAR on your data... One small issue here is the multicores, i should multithread but the rest of the script cannot be threaded... If you use nextflow this would have worked but i dont want to add another tool 
    """
    # Create the logger
    process = psutil.Process(os.getpid())
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('%s\tStarting Ppsyfinder, including The Alignment', time.ctime())
    CommandLine  = "%s --Method Fastq --pseudoCandidateDepth %s --InsertDistance %s --ChimericPairDepthTreshold %s --ChimericPairBinningTreshold %s --ChimericReadDepthTreshold %s --ChimericReadBinningTreshold %s --MergeChimReadWithChimpairTresh %s -STARindex %s -S %s -R1 %s -R2 %s" %(sys.argv[0], arguments.Psdepth,arguments.insdistance, arguments.ChimPairDepthTresh, arguments.ChimPairBinningTresh, arguments.ChimReadDepthTresh, arguments.ChimReadBinningTresh, arguments.chimreadpairdistance, arguments.STARindex,arguments.Sample,arguments.Fastq1, arguments.Fastq2)
    CheckingPaths()
    (genecoords,pseudogenecoords,exoncoords)=database()
    (Cleaninglist, MovingList,Outputfolder)=CreatingOutputDir(Sample)
    baminput=AlignmentWithSTAR(Fastq1,Fastq2,STARindex,Sample,Outputfolder)
    chrominfo=CreateTheChromosomeInfoFile(baminput,Sample,Cleaninglist)
    (Pseudogenecandidatesbed,Cigarbam)=Pseuodogenecandidates(Sample,baminput,exoncoords,pseudogenecoords,chrominfo,genecoords,Cleaninglist, MovingList,Outputfolder,Psdepth)
    bamChimericPairs=extractingchimericReads(Sample, baminput, Cleaninglist, MovingList,insdistance)
    PseudogeneCandidateChimbed=ChimericReadsOverlapwithPseudogeneCandidates(Sample, bamChimericPairs,Pseudogenecandidatesbed, Cleaninglist, MovingList,ChimPairDepthTresh,ChimPairBinningTresh)
    bamClipped=extractingClipped(Sample,baminput, Cleaninglist,MovingList) 
    clippedbed=ClippedBamtoCoord(bamClipped,Cleaninglist, MovingList)
    clippedmapping=clippedreadbinning(clippedbed,Cleaninglist, MovingList,ChimReadDepthTresh,ChimReadBinningTresh)
    clippedwithpseudogeneoverlap=intersectClippedPseudogeneCandidates(Sample, clippedmapping, Pseudogenecandidatesbed, Cleaninglist, MovingList)
    (clipandchimevidence, DetectedPseudogenesList)=CombiningClippWithChimericReads(Sample, clippedwithpseudogeneoverlap, PseudogeneCandidateChimbed, Cleaninglist, MovingList,chimreadpairdistance)
    CountKnownPseudogenes(Sample,baminput,pseudogenecoords, MovingList)
    summaryOutAnnotatedBoth=AnnotateFusionPoint(Sample,clipandchimevidence,exoncoords,genecoords,Cleaninglist,MovingList)
    GvizPlottingForOutput(Sample,baminput,Cigarbam,summaryOutAnnotatedBoth,exoncoords,Cleaninglist,MovingList)
    cleaning(Cleaninglist, MovingList,Outputfolder)
    logging.info('%s\tPpsy Finished, results in output folder',time.ctime().split(" ")[-2])
    logging.info('Memory Consumed in Bytes:\t%s',(process.memory_info()[0]))
    logging.info('RunTime in Sec:\t%s', int(time.time()-start))



if __name__=='__main__':
  arguments=parseArgs()
  if arguments.method == "Bam":
      mainBam(arguments.baminput, arguments.Sample,arguments.Psdepth, arguments.insdistance,arguments.ChimPairDepthTresh,arguments.ChimPairBinningTresh,arguments.ChimReadDepthTresh,arguments.ChimReadBinningTresh, arguments.chimreadpairdistance)

  else: # The input is the Fastq method so you run star within the pipeline 
      mainFastq(arguments.Fastq1,arguments.Fastq2,arguments.STARindex, arguments.Sample,arguments.Psdepth, arguments.insdistance,arguments.ChimPairDepthTresh,arguments.ChimPairBinningTresh,arguments.ChimReadDepthTresh,arguments.ChimReadBinningTresh, arguments.chimreadpairdistance)
