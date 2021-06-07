#!/usr/bin/python
# -*- coding: utf-8 -*-



import sys
import os
from os import listdir
import argparse
import xlsxwriter


def parseArgs():
    parser = argparse.ArgumentParser(description='Creates the Excel summary report')
    parser.add_argument('-I','--input', dest='SampleFolder',help='Name of the output folder from Ppsy that you want to create the IGV for (required)',required=True)
    parser.add_argument('-O','--OutputFolder', dest='outfolder', help='Output folder that will contain the IGV plots (required)', required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def writeOutputFolder(outfolder): 
    """
    This part writes the output folder 
    """
    commandmakeout = "mkdir -p %s" %outfolder
    os.system(commandmakeout)




def writeIGVscript(outfolder,SampleFolder):
    """ 
    This reads the PPsy out report to write the IGV image 
    """

    Sample = SampleFolder.split("_PPsyOut")[0].split("/")[-1]
    Clippedbam=os.path.abspath(SampleFolder)+"/Alignments/" + Sample + ".Clipped.bam"
    Pairbam=os.path.abspath(SampleFolder)+"/Alignments/" + Sample + ".ChimericReadPairs.bam"
   # print(os.path.abspath(SampleFolder))
    #print(os.listdir(SampleFolder+"/PpsyReports/"))
    report=os.path.abspath(SampleFolder)+"/PpsyReports/"+ Sample +".ChimPairs_ChimReads.Ppsy.txt"
    counter=0
    with open(report, "r") as infile: 
        next(infile)
        for l in infile:
            counter+=1
            l=l.strip()
            chimreads=l.split("\t")[13]
            chimpairs=l.split("\t")[5]
            chrompair=l.split("\t")[4]
            chromreads=l.split("\t")[12]
            pseudogene=l.split("\t")[0]
            if not chimreads == "NA": # If we have the reads use them as the breakpoint 
                rangestart=int(chimreads)-int(500)
                rangeend=int(chimreads)+500
                outfile=outfolder+"/"+pseudogene+"in"+str(chromreads)+"_"+str(chimreads)+"_"+str(counter)+".txt"
                outfilepng=outfolder+"/"+pseudogene+"in"+str(chromreads)+"_"+str(chimreads)+"_"+str(counter)+".png"
                with open(outfile, "w") as out: 
                    print("new\ngenome hg19", file=out)
                    print("load ",Pairbam, file=out )
                    print("load ",Clippedbam, file=out )
                    print("goto %s:%s-%s" %(str(chromreads),str(rangestart),str(rangeend)), file=out)
                    print("snapshot ",outfilepng,file=out)
                    print("exit", file=out)
            else: # Then whe use the start of chimpairleft
                rangestart=int(chimpairs)-int(500)
                rangeend=int(chimpairs)+500
                outfile=outfolder+"/"+pseudogene+"in"+str(chrompair)+"_"+str(chimpairs)+"_"+str(counter)+".txt"
                outfilepng=outfolder+"/"+pseudogene+"in"+str(chrompair)+"_"+str(chimpairs)+"_"+str(counter)+".png"
                with open(outfile, "w") as out: 
                    print("new\ngenome hg19", file=out)
                    print("load ",Pairbam, file=out )
                    print("load ",Clippedbam, file=out )
                    print("goto %s:%s-%s" %(str(chrompair),str(rangestart),str(rangeend)), file=out)
                    print("snapshot ",outfilepng,file=out)
                    print("exit",file=out)

def plotIGV(outfolder):
    """
    This script plots the files igv files in the output folder 
    """
    outf=os.path.abspath(outfolder)
    for i in os.listdir(outf):
        command="igv.sh --batch=" + outf +"/"+i 
        print(command)
        os.system(command)
        

def main(SampleFolder,outfolder):
    writeOutputFolder(outfolder)
    print("Writing script for IVG ...")
    writeIGVscript(outfolder,SampleFolder)
    #WriteTableWithoutLinks(outfolder,SampleFolders)
    print("plotting IGV")
    plotIGV(outfolder)
    print("Done!")


if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.SampleFolder,arguments.outfolder)
