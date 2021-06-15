import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()


setuptools.setup(
	name='DATtoRADMC',
	version='1.3',
	scripts=['DATtoRADMC'] ,
    author="Pablo Rutschmann",
    author_email="rupablo@ethz.ch",
    description="A .dat to RADMC3d file conversion class",
    long_description=long_description,
   	long_description_content_type="text/markdown",
    url="https://github.com/pablorutschmann/DATtoRADMC",
    packages=setuptools.find_packages(),
    classifiers=[
         "Programming Language :: Python :: 3.8.8",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
         ],
	