#!/usr/bin/python
# -*- coding: utf-8 -*-
# Idea is that the user input is all the folders containing the output. For know the folder itself is called the sample name. The putput is a folder named PpsySummary With the htlm within 
# To think About!  
# You need to pinpoint towards the folder _outputPseudogeneSearch! Before the _ is always the sample name! 


import sys
import os
from os import listdir
#from os.path import isfile, join

def readingFiles():
    SampleFolders=sys.argv[1:]
    # Create output Folder 
    commandmakeout = "mkdir -p PpsySummary"
    os.system(commandmakeout)
    return SampleFolders

def Writingtheccs():
    with open("PpsySummary/a.css", "w") as outccs:
        print >> outccs, """

.withBorders{
border-collapse: collapse;
 width: 100%;
 }

.Columnsnames{
font-size:150%;
text-align: center;
}

p{
font-family: Verdana,Sans-serif;
}

h1{
font-family: Verdana,Sans-serif;
}

table, th, td {
    border: 1px solid black;
    text-align: center;
    font-family: Verdana,Sans-serif; 
}

tr:nth-child(even) {
    background-color: #dce5ef;
}

        """
        

def WritingMainTable():
    with open("PpsySummary/PpsySummaryReport.html", "w") as outhtml:
        print >> outhtml, """
<!DOCTYPE html>
<script src="sorttable.js"></script>
<html>
<head>
<link rel="stylesheet" href="a.css">
<title>P&Psi;Finder Out</title>
</head>
<body>
<h1><p>P&Psi;Finder Summary Report</p></h1>
<hr> 
<p>Summary from P&Psi;Finder</p>
<table id="PseudoTable" class="withBorders sortable">
    <tr class = "Columnsnames"> 
        <th>Sample</th>
        <th>Report</th>
        <th><p>Known P&Psi;</p></th>
        <th><p>P&Psi;</p></th>
        <th>ChimRead Count</th>
        <th>ChimPair Count</th>
        <th>Insert Site</th>
        <th>Plot</th>
    </tr>

        """
def WritingTheRows(SampleFolders):
    for i in SampleFolders: 
        ReportFolder = i + "/PpsyReports"
        PlottingFolder = i + "/Plotting"
        #CircosFolder = i + "/Circos/CircosPictures"
        #CircosFolder = os.path.abspath(CircosFolder)
        with open("PpsySummary/PpsySummaryReport.html","a") as outhtml:
            # Now list all the files in the outputdir of ppsyFinder 
            filesinReportdir = os.listdir(ReportFolder)
            #filesinCircosdir = os.listdir(CircosFolder)
            filesinPlotdir = os.listdir(PlottingFolder)
            Exactfile = ReportFolder + "/" + "".join([s for s in filesinReportdir if ".ChimPairs_ChimReads.Ppsy.txt" in s]) # Get the exact coord file and save the string to variable open this file and loops thorugh it, first loop print link to known pseuodegenes and outlog, if the other loop skipt to print this
            KnownPpsy = "../" + ReportFolder + "/" + "".join([s for s in filesinReportdir if ".KnownProcessedPseudogenes.Out" in s])
            Sample = i.split("_PPsyOut")[0].split("/")[-1]
            
            ### This will be changed! Use a user defined outputName for all your Files 
            SampleNameForAlignmentFile =  Sample + ".r1.fq.gzAligned.sortedByCoord.out" #Obs change This So All intermediate Files will use the user named OutPut file!! 
            
            if os.path.isfile(Exactfile): # If the exact coords are there, Write this to the table!
                with open(Exactfile, "r") as f: # This part return the amount of rows in the file containing the pseudogenes 
                    for i, l in enumerate(f):
                        pass
                    if i == 1: # obs the file only contains the header, so it is empty then you break it here 
                        print >> outhtml, """
                    <tr>
           <td>%s</td>
           <td>-</td>
           <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p></td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
         </tr>
                    """ %(Sample, KnownPpsy)
                        
                        break
                    else:
                        amountofrows = int(i) + 1             
                with open(Exactfile, "r") as exactcoords:
                    next(exactcoords)# Skip the first row as it contains the header
                    rowintable=0 # Counter to make sure that it is only the first row in the sample that contains the KnownpseudogeneReport and the output Report
                    linecounter = 0 # this counter works as a control point, if you dont have any pseudogenes that dont contains evidence at both the chimpairs and chimreads we need to break the loop.
                    for line in exactcoords: # OBS!!! This is very dependent on the exact fusion coord! if you change something there you need to keep track and change it here aswell!
                        linecounter += 1
                        strippedline = line.strip()
                        Pseudogene = strippedline.split("\t")[0]
                        Exactfile = "../" + ReportFolder + "/" + "".join([s for s in filesinReportdir if ".ChimPairs_ChimReads.Ppsy.txt" in s]) # When you create it the second time you need to write as it is inserted into the html report 
                        insertsite = strippedline.split("\t")[12] + ":" + str(strippedline.split("\t")[13]) + "-" + str(strippedline.split("\t")[14])
                        # circospicture = strippedline.split("\t")[0] + "_" + strippedline.split("\t")[1] + "_" + strippedline.split("\t")[2] + ".png"
                        chimreadCounts = str(strippedline.split("\t")[7])
                        chimpairCounts = str(strippedline.split("\t")[15])
                        #circospictureend = "%s_%s_%s_%s_%s_%s_%s_%s_%s.png" %(Pseudogene,strippedline.split("\t")[1],  str(strippedline.split("\t")[2]), str(strippedline.split("\t")[3]),  str(strippedline.split("\t")[4]),  str(strippedline.split("\t")[5]),  str(strippedline.split("\t")[6]),str(strippedline.split("\t")[7]),str(strippedline.split("\t")[8]),str(strippedline.split("\t")[9]),str(strippedline.split("\t")[10]),str(strippedline.split("\t")[11]),str(strippedline.split("\t")[12]),str(strippedline.split("\t")[13]),str(strippedline.split("\t")[14]),str(strippedline.split("\t")[15]))
                        if not strippedline.split("\t")[12] == "NA": 
                            PlottingPictureend = "%swithin%s-%s.png" %(Pseudogene,strippedline.split("\t")[12],strippedline.split("\t")[13])
                        else: 
                            PlottingPictureend = "%swithin%s-%s.png" %(Pseudogene,strippedline.split("\t")[4],strippedline.split("\t")[5])

                        print PlottingPictureend
                        #circospicture  = "".join([s for s in filesinCircosdir if circospictureend in s])

                        PlottingPicture="".join([s for s in filesinPlotdir if PlottingPictureend in s])
                        #CircosPicturePath = "../" + CircosFolder + "/" + circospicture # One step back from teh ppsy report, then you can use the path
                        PlottingPicturePath =  "../" + PlottingFolder + "/" + PlottingPicture 
                        print PlottingPicturePath
                        
                        #print CircosFolder
                        insertanno = strippedline.split("\t")[16]
                        # Here we change the coloring depending on where the insertsite is, if intergenic green, if intronic blue, if exonic red 
                        if insertanno == "exonic":
                            anncolor = "#990000"
                        
                        elif insertanno == "intronic":
                            anncolor = "#0066ff"
                        
                        else: 
                            anncolor = "#006600"
                        if not chimreadCounts == "NA" and not chimpairCounts == "NA":
                            rowintable+=1
                            if rowintable == 1:
                                print >> outhtml, """
    <tr>
           <td>%s</td>
           <td><a href=\"%s\", target="_blank">Output</a></td>
           <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p> </td>
           <td><i>%s<i></td>
           <th>%i</th>
           <th>%i</th>
           <td style=\"color:%s;\"><p>&rarr; %s</p></td>       
           <td><a href=\"%s\", target=\"_blank\"> <img src=\"%s\" style=\"width:42px;height:42px;border:0\"></a> </td>
         </tr>
                            """ %(Sample, Exactfile, KnownPpsy, Pseudogene, int(chimreadCounts),int(chimpairCounts),anncolor,insertsite,PlottingPicturePath,PlottingPicturePath)

                            else: # If we have done the first line for the Psedodogene skip writing the first three columns 
                                print >> outhtml, """
    <tr>
           <td>%s</td>
           <td><a href=\"%s\", target="_blank">Output</a></td>
           <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p></td>
           <td><i>%s</i></td>
           <th>%i</th>
           <th>%i</th>
           <td style=\"color:%s;\"><p>&rarr; %s</p></td>
           <td><a href=\"%s\", target=\"_blank\"> <img src=\"%s\" style=\"width:42px;height:42px;border:0\"></a> </td>
         </tr>
                            """ %(Sample,Exactfile,KnownPpsy,Pseudogene, int(chimreadCounts),int(chimpairCounts),anncolor,insertsite, PlottingPicturePath,PlottingPicturePath)

                        elif rowintable < 1 and linecounter+1 == amountofrows: # This is because we dont hit any pseudogene with evidence from both. Therefore we wont print them in the report
                            print >> outhtml, """
                    <tr>
           <td>%s</td>
           <td><a href=\"%s\", target="_blank">Output</a></td>
           <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p></td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
           <td>-</td>
         </tr>

                    """ %(Sample,Exactfile, KnownPpsy)                            
                            
            else: # Samples does not have any pseudogenes  
                print >> outhtml, """
                <tr>
       <td>%s</td>
       <td><a href=\"%s\", target="_blank">Output</a></td>         
       <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p></td>
       <td>-</td>
       <td>-</td>
       <td>-</td>
       <td>-</td>
       <td>-</td>
     </tr>
                
                """ %(Sample,Exactfile,KnownPpsy)
        
    with open("PpsySummary/PpsySummaryReport.html","a") as outhtml:
        print >> outhtml, """
        </table>
   
        
        </body>
        </html>
        """

def main():
    print "Parsing the html report ..."
    SampleFolders=readingFiles()
    Writingtheccs()
    WritingMainTable()
    WritingTheRows(SampleFolders)
    print "Done!"
main()
