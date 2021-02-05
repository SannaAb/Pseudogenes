#!/usr/bin/python
# -*- coding: utf-8 -*-


import sys
import os
from os import listdir
import datetime
import glob
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='Creates the output html from the PPSY run, you descide what you want to display, pseudogenes detected by both evidence chimpairs and chimreads or is it enough with only evidence from one?')
    parser.add_argument('-I','--inputs', dest='SampleFolders',nargs='+',help='Name of the output folder/folders from Ppsy that you want to create the report for (required)',required=True)
    parser.add_argument('-O','--OutputFolder', dest='outfolder', help='Output folder that contains the html report (required)', required=True)
    parser.add_argument("-f", "--format", type=str, dest='form',choices=['strict','loose'],help='What will the output report contain? If you choose strict only output hits that are supported by both chimreads and chim pairs are presented if you choose loose all the pseudoenge candidates with their insertsites are presented (required)', default='loose')
    arguments=parser.parse_args(sys.argv[1:])
    return arguments

def readingFiles():
    SampleFolders=sys.argv[1:]
    # Create output Folder 
    commandmakeout = "mkdir -p PpsySummary"
    os.system(commandmakeout)
    return SampleFolders

def Writingtheccs(outfolder):
    commandmakeout = "mkdir -p %s" %outfolder
    os.system(commandmakeout)
    css = outfolder + "/a.css"
    with open(css, "w") as outccs:
        print("""

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
    #border: 1px solid black;
    border:none;
    text-align: center;
    font-family: Verdana,Sans-serif; 
}

tr:nth-child(even) {
    background-color: #dce5ef;
}

.tooltip {
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black;
}

.tooltip .tooltiptext {
  visibility: hidden;
  width: 250px;
  background-color: black;
  color: #fff;
  text-align: center;
  border-radius: 6px;
  padding: 5px 0;
    
  /* Position the tooltip */
  position: absolute;
  z-index: 1;
  bottom: 100%;
  left: 50%;
  margin-left: -125px;
}



.tooltip:hover .tooltiptext {
  visibility: visible;
  font-size: 10px;
}



/* basic positioning */
.legend { list-style: none; }
.legend li { float: left; margin-right: 10px; }
.legend span { border: 1px solid #ccc; float: left; width: 12px; height: 12px; margin: 2px; }
/* your colors */
.legend .intergenic { background-color: #000000; }
.legend .intronic { background-color: #198b19; }
.legend .exonic { background-color: #990000; }

        """, file=outccs)
        

def WritingMainTable(outfolder,SampleFolders):
    ppsyhtml = outfolder + "/PpsySummaryReport.html"
    date = str(datetime.datetime.now()).split(":")[0]+":"+str(datetime.datetime.now()).split(":")[1]
    nsamples = len(SampleFolders)
    nsampleswithpseudogene = 0 
    path = "/PpsyReports/*ChimPairs_ChimReads.Ppsy.txt"
    fullpathtoppsyreport = [s + path for s in SampleFolders] # here we append the full path, just to see how many pseudogenes we detect

    for i in fullpathtoppsyreport:
        textfile = glob.glob(i)
        textfile= "".join(textfile)
        try: 
            with open(textfile, "r") as report: 
                if sum(1 for _ in report) > 1: 
                    nsampleswithpseudogene +=1 
        except IOError: 
            # We dont have the report
            continue
    with open(ppsyhtml, "w") as outhtml:
        print("""
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
<p>Time for report generation: <strong>%s</strong>
<br>Screening nSamples: <strong>%s</strong>
<br>Detected Pseudogenes in nSamples: <strong>%s</strong> 
<br>
</p>

<p>Insert annotation color key:
<ul class="legend">
    <li><span class="intergenic"></span> Intergenic</li>
    <li><span class="intronic"></span> Intronic</li>
    <li><span class="exonic"></span> Exonic</li>
</ul>
</p>
<table id="PseudoTable" class="withBorders sortable">
    <tr class = "Columnsnames"> 
        <th>
        <div class="tooltip">Sample
  <span class="tooltiptext">Analyzed sample</span>
        </div></th>
        
        <th>
        <div class="tooltip">Report
  <span class="tooltiptext">Ppsy report for the sample</span>
        </div></th>
        <th><p><div class="tooltip"> Known P&Psi;<span class="tooltiptext">Known Ppsy report for the sample</span></div></p> </th>
        <th><p><div class="tooltip"> P&Psi;<span class="tooltiptext">Detected Pseudogene for the sample</span></div><p></th>
        <th><div class="tooltip">ChimPair Count<span class="tooltiptext">Amount of left chimeric pairs supporting the insertsite for the specific pseudogene</span></div></th>
        <th><div class="tooltip">ChimRead Count<span class="tooltiptext">Amount of chimeric reads supporting the insertsite for the specific pseudogene</span></div></th>
        <th><div class="tooltip">Insert Site<span class="tooltiptext">Insertsite of the Pseudogene, if chimeric reads are not available the insert point is from left chimeric pairs instead</span></div></th>
        <th><div class="tooltip">Plot<span class="tooltiptext" style="margin-left: -225px;">Coverage plot of the Pseudogene and its insert</span></th>
    </tr>

        """ % (date,nsamples,nsampleswithpseudogene), file=outhtml)

#<tr><td style="border-top:1px solid black;" colspan="8"></td></tr>

def WritingTheRows(SampleFolders,outfolder,form):
    print("Parsing the html Report ...")
    for m in SampleFolders: 
        #print "Parsing\t%s" %m
        ReportFolder = m + "/PpsyReports"
        PlottingFolder = m + "/Plotting"
        outhtmlfile = outfolder + "/PpsySummaryReport.html"
        with open(outhtmlfile,"a") as outhtml:
            # Now list all the files in the outputdir of ppsyFinder 
            filesinReportdir = os.listdir(ReportFolder)
            #filesinCircosdir = os.listdir(CircosFolder)
            filesinPlotdir = os.listdir(PlottingFolder)
            Exactfile = ReportFolder + "/" + "".join([s for s in filesinReportdir if ".ChimPairs_ChimReads.Ppsy.txt" in s]) # Get the exact coord file and save the string to variable open this file and loops thorugh it, first loop print link to known pseuodegenes and outlog, if the other loop skipt to print this
            KnownPpsy = "../" + ReportFolder + "/" + "".join([s for s in filesinReportdir if ".KnownProcessedPseudogenes.Out" in s])
            Sample = m.split("_PPsyOut")[0].split("/")[-1]
            
            ### This will be changed! Use a user defined outputName for all your Files 
            SampleNameForAlignmentFile =  Sample + ".r1.fq.gzAligned.sortedByCoord.out" #Obs change This So All intermediate Files will use the user named OutPut file!! 
            i=0 # Count amount of rows 
            if os.path.isfile(Exactfile): # If the exact coords are there, Write this to the table!
                
                
                with open(Exactfile, "r") as f: # This part return the amount of rows in the file containing the pseudogenes 
                    for i, l in enumerate(f):
                        pass
                    amountofrows=int(i)+1                    
                    if amountofrows == 1: # obs the file only contains the header, so it is empty then you break it here 
                        print("""
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
                        """ %(Sample, KnownPpsy), file=outhtml)
                        
                        continue
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
                        insertsitechimpair = strippedline.split("\t")[4] + ":" + str(strippedline.split("\t")[5]) + "-" + str(strippedline.split("\t")[6])
                        chimreadCounts = str(strippedline.split("\t")[7])
                        chimpairCounts = str(strippedline.split("\t")[15])
                        if not strippedline.split("\t")[12] == "NA": 
                            PlottingPictureend = "%swithin%s-%s.png" %(Pseudogene,strippedline.split("\t")[12],strippedline.split("\t")[13])
                        else: 
                            PlottingPictureend = "%swithin%s-%s.png" %(Pseudogene,strippedline.split("\t")[4],strippedline.split("\t")[5])
                        PlottingPicture="".join([s for s in filesinPlotdir if PlottingPictureend in s])
                        PlottingPicturePath =  "../" + PlottingFolder + "/" + PlottingPicture 
                        insertanno = strippedline.split("\t")[16]
                        # Here we change the coloring depending on where the insertsite is, if intergenic black, if intronic green, if exonic red 
                        if insertanno == "exonic":
                            anncolor = "#990000"
                        elif insertanno == "intronic":
                            anncolor = "#198b19"
                        else: 
                            anncolor = "#000000"
                        if form == 'strict':
                            if not chimreadCounts == "NA" and not chimpairCounts == "NA":
                                rowintable+=1
                                if rowintable == 1:
                                    print("""
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
                                    """ %(Sample, Exactfile, KnownPpsy, Pseudogene, int(chimreadCounts),int(chimpairCounts),anncolor,insertsite,PlottingPicturePath,PlottingPicturePath), file=outhtml)

                                else: # If we have done the first line for the Psedodogene skip writing the first three columns 
                                    print("""
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
                                    """ %(Sample,Exactfile,KnownPpsy,Pseudogene, int(chimreadCounts),int(chimpairCounts),anncolor,insertsite, PlottingPicturePath,PlottingPicturePath), file=outhtml)

                            elif rowintable < 1 and linecounter+1 == amountofrows: # This is because we dont hit any pseudogene with evidence from both. Therefore we wont print them in the report
                                print("""
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
                                """ %(Sample,Exactfile, KnownPpsy), file=outhtml)                            
                        elif form == 'loose': # We chose to print everything to the ppsy report
                            if insertsite == "NA:NA-NA": # The clipped reads are NA by evidence so we go for chimpair left anchor instead as the insert coords in the report
                                insertsite=insertsitechimpair
                            rowintable+=1
                            if rowintable == 1:
                                print("""
        <tr>
               <td>%s</td>
               <td><a href=\"%s\", target="_blank">Output</a></td>
               <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p> </td>
               <td><i>%s<i></td>
               <th>%s</th>
               <th>%s</th>
               <td style=\"color:%s;\"><p>&rarr; %s</p></td>       
               <td><a href=\"%s\", target=\"_blank\"> <img src=\"%s\" style=\"width:42px;height:42px;border:0\"></a> </td>
             </tr>
                                """ %(Sample, Exactfile, KnownPpsy, Pseudogene, chimreadCounts,chimpairCounts,anncolor,insertsite,PlottingPicturePath,PlottingPicturePath), file=outhtml)

                            else: # If we have done the first line for the Psedodogene skip writing the first three columns 
                                print("""
        <tr>
               <td>%s</td>
               <td><a href=\"%s\", target="_blank">Output</a></td>
               <td><p><a href=\"%s\", target="_blank">KnownP&Psi;</a></p></td>
               <td><i>%s</i></td>
               <th>%s</th>
               <th>%s</th>
               <td style=\"color:%s;\"><p>&rarr; %s</p></td>
               <td><a href=\"%s\", target=\"_blank\"> <img src=\"%s\" style=\"width:42px;height:42px;border:0\"></a> </td>
             </tr>
                                """ %(Sample,Exactfile,KnownPpsy,Pseudogene, chimreadCounts,chimpairCounts,anncolor,insertsite, PlottingPicturePath,PlottingPicturePath), file=outhtml)

                        elif rowintable < 1 and linecounter+1 == amountofrows: # This is because we dont hit any pseudogene with evidence from both. Therefore we wont print them in the report
                            print("""
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
                            """ %(Sample,Exactfile, KnownPpsy), file=outhtml)


            else: # Samples does not have any pseudogenes  
                print("""
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
                
                """ %(Sample,Exactfile,KnownPpsy), file=outhtml)
        
    with open(outhtmlfile,"a") as outhtml:
        print("""
        </table>
   
        
        </body>
        </html>
        """,file=outhtml)

def main(SampleFolders,outfolder,form):
    Writingtheccs(outfolder)
    WritingMainTable(outfolder,SampleFolders)
    WritingTheRows(SampleFolders,outfolder, form)
    print("Done!")


if __name__=='__main__':
  arguments=parseArgs()
  main(arguments.SampleFolders,arguments.outfolder,arguments.form)


