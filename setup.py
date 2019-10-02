import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
     name='Ppsy',  
     version='0.1',
     scripts=['Ppsy.py'] ,
     author="Sanna Abrahamsson",
     author_email="sannaabrahamsson@gmail.com",
     description="Script to detect Processed Pseudogenes",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/SannaAb/Pseudogenes",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 2", 
         
     ],
 )
