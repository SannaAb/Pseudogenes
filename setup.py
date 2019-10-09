#!/usr/bin/py 
import setuptools


setuptools.setup(
    name='Ppsy',  
    version='0.1',
    scripts=['Scripts/Ppsy.py','Scripts/MakePPsyReport.py'] ,
    author="Sanna Abrahamsson",
    author_email="sannaabrahamsson@gmail.com",
    description="Script to detect Processed Pseudogenes",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SannaAb/Pseudogenes",
    packages=setuptools.find_packages(),
    
    package_data={'Ppsy': ['README.md', 'HG19_databases/Exon_coord_hg19_refgene.bed',"HG19_databases/Gene_coord_hg19_refgene.bed","HG19_databases/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed","HG19_databases/Ppsy_Config.txt"]
                   },
    classifiers=[
         "Programming Language :: Python :: 2", 
         
     ],
 )
