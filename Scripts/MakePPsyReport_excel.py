#!/usr/bin/python
# -*- coding: utf-8 -*-



import sys
import os
from os import listdir
#from os.path import isfile, join
import argparse
import xlsxwriter
#pip install XlsxWriter


def parseArgs():
    parser = argparse.ArgumentParser(description='Creates the Xcell summary reports, one with out links and one with')
    parser.add_argument('-I','--inputs', dest='SampleFolders',nargs='+',help='Name of the output folder/folders from Ppsy that you want to create the report for (required)',required=True)
    parser.add_argument('-O','--OutputFolder', dest='outfolder', help='Output folder that will contain the Excell reports (required)', required=True)
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def writeOutputFolder(outfolder): 
    """
    This part writes the output folder 
    """
    commandmakeout = "mkdir -p %s" %outfolder
    os.system(commandmakeout)

def WriteTableWithoutLinks(outfolder,SampleFolders): 
    """
    This report writes without links so it can be easy passable, no need to keep the path structure
    """
    outfile = outfolder + "/Ppsy_summary_withoutLinks.xlsx"
    woorkbook = xlsxwriter.Workbook(outfile)
    worksheet = woorkbook.add_worksheet()
    rownumber = 0
    bold = woorkbook.add_format({'bold': True})
    # Set the column names
    columns = ('Sample','Pseudogene','ChimPair Count','ChimRead Count','Insert Site', 'region', 'annotation')
    worksheet.write_row(rownumber,0, columns, bold)
    rownumber = 0 # We start at this row
    for m in SampleFolders: # Lets loop the sample folders, everything for this report will be in the report of the excell files
        Sample = m.split("_PPsyOut")[0].split("/")[-1]
        print "Parsing\t%s" %Sample
        ReportFolder = m + "/PpsyReports"
        filesinReportdir = os.listdir(ReportFolder)
        Exactfile = ReportFolder + "/" + "".join([s for s in filesinReportdir if ".ChimPairs_ChimReads.Ppsy.txt" in s]) # Get the exact coord file and save the string to variable open this file and loops thorugh it, first loop print link to known pseuodegenes and outlog, if the other loop skipt to print this
        if os.path.isfile(Exactfile): # If the exact coords are there, Write this to the table ,it should always be though.. Otherwise it should print a file with just the header
                with open(Exactfile, "r") as f: # This part return the amount of rows in the file containing the pseudogenes 
                    for i, l in enumerate(f):
                        pass
                    amountofrows=int(i)+1
                    if amountofrows == 1: # obs the file only contains the header, so it is empty then you break it here 
                        columns = (Sample,'-','-','-','-','-','-')
                        rownumber+=1
                        worksheet.write_row(rownumber,0, columns)
                    else: # We have detected pseudogenes in the others, write them to the excell 
                        f.seek(0) # reset the loop
                        next(f) # skip the first column
                        for l in f: 
                            rownumber+=1
                            l=l.strip()
                            Pseudogene = l.split("\t")[0]
                            leftchimpaircount = l.split("\t")[7]
                            rightchimpaircount = l.split("\t")[11]
                            if not leftchimpaircount == "NA" and not rightchimpaircount == "NA": # if none of the pairs are NA add them togehter
                                ChimPairCount = int(leftchimpaircount) + int(rightchimpaircount) 
                            elif not leftchimpaircount == "NA" and rightchimpaircount == "NA": # else only use left pair 
                                ChimPairCount = int(leftchimpaircount)
                            else: # if left is NA set it as NA
                                ChimPairCount = "NA"
                            ChimReadCount = l.split("\t")[15]
                            if not ChimReadCount == "NA": 
                                ChimReadCount = int(ChimReadCount)
                            # Insertsite, if we have the chim read split use that, otherwise use the left chimpair count
                            if l.split("\t")[12] != "NA":
                                insertsite = l.split("\t")[12] + ":" + str(l.split("\t")[13]) + "-" + str(l.split("\t")[14])
                            else: 
                                insertsite = l.split("\t")[4] + ":" + str(l.split("\t")[5]) + "-" + str(l.split("\t")[6])
                            insertregion = l.split("\t")[16]
                            insertranno = l.split("\t")[-1]
                            columns = (Sample,Pseudogene,ChimPairCount,ChimReadCount,insertsite, insertregion,insertranno)
                            worksheet.write_row(rownumber,0, columns)
        else: 
            print Sample
    woorkbook.close()

    

def main(SampleFolders,outfolder):
    writeOutputFolder(outfolder)
    print "Parsing the xcell report without hyperlinks ..."
    WriteTableWithoutLinks(outfolder,SampleFolders)
    print "Done!"


if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.SampleFolders,arguments.outfolder)
