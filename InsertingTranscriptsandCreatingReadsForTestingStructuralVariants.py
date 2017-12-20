#!/usr/bin/python

import sys
import subprocess
import os
import argparse
from Bio import SeqIO 

# module load bedtools
# module load emboss

def parseArgs(argv):
    '''
    parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Extracts all exons of the desired gene from the reference genome, merge overlapping coordinates. Insert that gene somewhere else in the genome. Then create reads in that area where the transcript is inserted')
    parser.add_argument("-g", dest = 'geneofinterest',help ="The gene that you want to insert")
    parser.add_argument("-i", dest = 'insertcoords',help ="The coordinates where you want to insert the transcript, ex.. chr1:10004")
    parser.add_argument("-dup", dest = 'exonduplication',help ="add which exon that should be duplicated, example exon 1, insert after exon 3 -dup 1:3", required=False)
    parser.add_argument("-d", dest = 'exondeletion',help ="add which exon that should be deleted, example exon 1 -d 1", required=False)
    parser.add_argument("-s", dest = 'swapbetweenexons',help ="add which exons you want to swap place between, example exon 1 and 2 -s 1:2", required=False)
    parser.add_argument("-r", dest = 'reversecomplement',help ="add which exon that should set to reverse, example exon 1 -r 1", required=False)
    parser.add_argument("-rent", dest = 'reverseentire',help ="ad this flag if you want to reverse the entire gene insert, set rent TRUE", required=False)
    arguments = parser.parse_args(argv)
    return arguments


def getTranscript(geneofinterest,exonduplication,exondeletion,swapbetweenexons,reversecomplement,reverseentire):
    bedfileexoncoords = geneofinterest + ".bed.tmp"
    fastafileexoncoords = geneofinterest + ".fasta.tmp"
    exoncoordsCommand = "grep -w \"%s\" /jumbo/WorkingDir/B17-006/db/Exongene_Homo_sapiens.GRCh37_Refseq.bed | awk -v OFS='\t' '{print $1,$2-1,$3,$4}' |sort -k1,1 -k2,2n | mergeBed > %s" % (geneofinterest,bedfileexoncoords) # Greps the exon coords of the wanted gene, merges overlapping exons
    os.system(exoncoordsCommand)
   # exoncoordsCapture = subprocess.Popen(exoncoordsCommand, shell=True, stdout=subprocess.PIPE) # Subprocess Pipe
   # exoncoordsForGeneOfInterest = exoncoordsCapture.stdout.read()
   # print exoncoordsForGeneOfInterest
    getSequenceFromFasta = "bedtools getfasta -fi /jumbo/db/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_nuHg19_mtrCRS.fa -bed %s -fo %s" % (bedfileexoncoords, fastafileexoncoords) # Extract the fasta sequences based on the coordinates from the merged exon bedfile
    #print getSequenceFromFasta
    os.system(getSequenceFromFasta)
    
    #---------------------------------------------------------------------------#
    
    # Exon Reconstruction
    
    #---------------------------------------------------------------------------#

    # Check if you want an exon duplication

    if exonduplication is not None:
        exonduplicationindex_exonofinterest = int(exonduplication.split(":")[0])-1 
        exonduplicationindex_insert = int(exonduplication.split(":")[1])
        with open(fastafileexoncoords, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            records.insert(exonduplicationindex_insert,records[exonduplicationindex_exonofinterest] ) # duplicate the chosen exon and insert it at user defined position
        with open(fastafileexoncoords, "w") as outputhandle:
            SeqIO.write(records, outputhandle, "fasta")


        
    # Check if you want a deletion 
    if exondeletion is not None:
        exondeletionindex = int(exondeletion)-1
        with open(fastafileexoncoords, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            del records[exondeletionindex]
        with open(fastafileexoncoords, "w") as outputhandle:
            SeqIO.write(records, outputhandle, "fasta")

    # Check if you want an swap between two exons
    
    if swapbetweenexons is not None:
        firstexonindex = int(swapbetweenexons.split(":")[0]) - 1 # Convert user input to python index in list, first item 
        secondexonindex = int(swapbetweenexons.split(":")[1])- 1 # Convert user input to python index in list, second item
        with open(fastafileexoncoords, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            records[firstexonindex],records[secondexonindex] = records[secondexonindex], records[firstexonindex] # Swapping the items in the records list
        with open(fastafileexoncoords, "w") as outputhandle:
            SeqIO.write(records, outputhandle, "fasta")   
    
    # Check if you want to reverse one of the exons

    if reversecomplement is not None:
        indexforexonreversecomplement = int(reversecomplement)-1
        with open(fastafileexoncoords, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            print records
            records[indexforexonreversecomplement] = records[indexforexonreversecomplement].reverse_complement()
            print records
          #  print records[indexforexonreversecomplement].reverse_comlement()
            #records[indexforexonreversecomplement] = records[indexforexonreversecomplement][::-1] # Reverse for the exon of interest
        with open(fastafileexoncoords, "w") as outputhandle:
            SeqIO.write(records, outputhandle, "fasta")

    
       
    CommandMergeexons = "grep -v \"^>\" %s | awk 'BEGIN { ORS=\"\"; print \">%s\\n\" } { print }'" %(fastafileexoncoords,geneofinterest) # Merging the fasta sequences to one big sequence 
    #print CommandMergeexons
    MergeexonsCapture = subprocess.Popen(CommandMergeexons, shell=True, stdout=subprocess.PIPE) # Subprocess Pipe
    MergeExonsForGeneOfInterest = MergeexonsCapture.stdout.read()
    #print MergeExonsForGeneOfInterest
    print "\n"

    # Check if you want to reverse the entire transcript
    if reverseentire is not None:
        with open("reverse.fasta.tmp", "w") as SequenceTorev:
            SequenceTorev.write(MergeExonsForGeneOfInterest)
        # Package revseq from emboss/6.6.0    
        commandReverseEntireGene = "revseq -sequence reverse.fasta.tmp -reverse -complement -outseq reverse2.fasta.tmp" # If you dont want the complement set nocomplement
        os.system(commandReverseEntireGene)
        with open("reverse2.fasta.tmp", "r") as reverseout:
            MergeExonsForGeneOfInterestlist = reverseout.readlines()
            MergeExonsForGeneOfInterest = MergeExonsForGeneOfInterestlist[0] + "".join(MergeExonsForGeneOfInterestlist[1:]).replace('\n', '') # removing new lines and pasting the first header of the fasta and the sequence back together
        #command_removetmp       
    # Remove the new bed and fasta file
    Command_remove = "rm %s %s" %(bedfileexoncoords,fastafileexoncoords)
    os.system(Command_remove)
    return MergeExonsForGeneOfInterest

def insertinreference(MergeExonsForGeneOfInterest, insertcoords, geneofinterest):
   # print MergeExonsForGeneOfInterest
    chromosomeforinsert = insertcoords.split(":")[0] # Tells us which chromosome we will insert our string 
    startPositionOfinsert = insertcoords.split(":")[1] # Tells us what the exact coordinate is for the insert transcript
    Commandgetsequenceforchromsomeofinterest = "samtools faidx /jumbo/db/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_nuHg19_mtrCRS.fa %s | grep -v \">\" | tr -d '\n'" % chromosomeforinsert
    FastaChromofinterestCapture = subprocess.Popen(Commandgetsequenceforchromsomeofinterest, shell = True, stdout=subprocess.PIPE)
    FastaChromofinterest = FastaChromofinterestCapture.stdout.read()
  #  print FastaChromofinterest

    
    newinsert = FastaChromofinterest[:int(startPositionOfinsert)] + str(MergeExonsForGeneOfInterest.split("\n")[1]) +FastaChromofinterest[int(startPositionOfinsert):] # Entire chromsome including the insert
    
    range500withinsertfasta = ">" + geneofinterest +  "_" + insertcoords +"\n" +  FastaChromofinterest[(int(startPositionOfinsert)-500):int(startPositionOfinsert)] + str(MergeExonsForGeneOfInterest.split("\n")[1]) + FastaChromofinterest[int(startPositionOfinsert):(int(startPositionOfinsert)+500)]
    

    # writing the new fasta, increasing range 500 to both sides so we can obtain the fusion point in the pseudogene script
    newinsertfastawith500bpindreasedrange = str(geneofinterest) + "_" + str(insertcoords) + "_insert.fasta"
    
    with open(newinsertfastawith500bpindreasedrange,"w") as f:
        f.write(range500withinsertfasta)  

    return newinsertfastawith500bpindreasedrange


def createreadpairs(newinsertfastawith500bpindreasedrange):
    basenameread = newinsertfastawith500bpindreasedrange.split(".")[0]
    replacebasenameread = basenameread.replace(":","_")# replaces : with _
    # Insert 
    Commanddwgsim = "wgsim -d 500 -N 100000 -R 0 -X 0 -r 0 %s -1 90 -2 90 %s.r1.fq %s.r2.fq" %(newinsertfastawith500bpindreasedrange, replacebasenameread,replacebasenameread) # Creates the read for the new sequence insert 
    os.system(Commanddwgsim)
    
def main(geneofinterest,insertcoords,exonduplication,exondeletion,swapbetweenexons,reversecomplement,reverseentire):

    
    MergeExonsForGeneOfInterest = getTranscript(geneofinterest,exonduplication,exondeletion,swapbetweenexons,reversecomplement,reverseentire)
    newinsertfastawith500bpindreasedrange=insertinreference(MergeExonsForGeneOfInterest, insertcoords,geneofinterest)
    createreadpairs(newinsertfastawith500bpindreasedrange)
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.geneofinterest, args.insertcoords,args.exonduplication, args.exondeletion, args.swapbetweenexons, args.reversecomplement, args.reverseentire)
