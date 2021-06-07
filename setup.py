#!/usr/bin/py 
import setuptools

setuptools.setup(
    name='Ppsy',  
    version='0.1',
    scripts=['Scripts/Ppsy.py','Scripts/MakePPsyReport_html.py', 'Scripts/MakePPsyReport_excel.py', 'Scripts/CreateCircosImageFromPpsyOut.py', 'Scripts/Ppsy.py', 'Scripts/MakePPsyReport_html.py', 'Scripts/CreateiGVImageFromPpsyOut.py'] ,
    author="Sanna Abrahamsson",
    author_email="sannaabrahamsson@gmail.com",
    description="Script to detect Processed Pseudogenes",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SannaAb/Pseudogenes",
    packages=setuptools.find_packages(),
    install_requires=['Cython','pandas','pysam','psutil','XlsxWriter'],
    package_data={'Ppsy': ['README.md', 'HG19_databases/Exon_coord_hg19_refgene.bed',"HG19_databases/Gene_coord_hg19_refgene.bed","HG19_databases/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed", "TestDataset/CBX3_chr15_40854180_insert.r1.fq","TestDataset/CBX3_chr15_40854180_insert.r2.fq"]
                   },
    data_files= [('HG19_databases',['HG19_databases/Exon_coord_hg19_refgene.bed','HG19_databases/Gene_coord_hg19_refgene.bed','HG19_databases/KnownProcessedPseudogenes_Homo_sapiens.GRCh37.75_CHR.bed', 'TestDataset/CBX3_chr15_40854180_insert.r1.fq','TestDataset/CBX3_chr15_40854180_insert.r2.fq'])],
    include_package_data=True,
    classifiers=[
         "Programming Language :: Python :: 3", 
         
     ],
 )
