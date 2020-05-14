DataGetter provides an interface for several MGEs detection softwares (Alien Hunter, IntegronFinder and MacSyFinder) that are based on different known MGEs features and complement one another.
It allows to run these three softwares on a dataset of genomes and returns the results as a tabulated file easily insertable into a database.

How to launch DataGetter :
./DataGetter.py <path to the directory containing the assemblies files> <path to the directory where to store the returned file>


How to install the three softwares :

############################################ MacSyFinder ############################################

wget https://bintray.com/gem-pasteur/MacSyFinder/download_file?file_path=macsyfinder-1.0.5.tar.gz

make sure dependencies are present
	- Python, (>=2.7 & <3.0)
	- makeblastdb (>= 2.2.28) or formatdb (>=2.2.26)
	- makeblastdb is a distribute with blastdb (http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
	- formatdb 
	- hmmsearch 3.1 (http://hmmer.janelia.org/)

tar xzvf download_file?file_path=macsyfinder-1.0.5.tar.gz
cd macsyfinder-1.0.5
python setup.py build
python setup.py
sudo python setup.py install

########################################### IntegronFinder ###########################################

make sure dependencies are present
	- Python 2.7
		- Pandas (0.18.0)
		- Numpy (1.9.1)
		- Biopython (1.65)
		- Matplotlib (1.4.3)
		- psutils (2.1.3)
	- HMMER (3.1b1)
	- INFERNAL (1.1)
	- Prodigal (2.6.2)

tar xzvf Integron_Finder-1.5.1.tar.gz
cd Integron_Finder-1.5.1
sudo pip install integron_finder

#add to .bashrc
export PATH=$PATH:/home/clem/Integron_Finder-1.5.1
# or # export INTEGRON_HOME=/home/clem/Integron_Finder-1.5.1

############################################ Alien-Hunter ############################################

make sure dependencies are present
	- Perl (5.6.1)
	- Java (1.4.2)

tar xzvf alien_hunter.tar.gz
cd alien_hunter-1.7
./alien_hunter
