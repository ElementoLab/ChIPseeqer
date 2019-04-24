==========================

ChIPseeqer (v2.0)

==========================


== Dependencies ==

ChIPseeqer has been tested in Linux Ubuntu distribution v10.04 LTS (Lucid Lynx), in x86/x86-64 
architectures.
	
1. To install the packages needed for Ubuntu, use apt-get. 

2. To download a package from apt-get, use the command:

sudo apt-get install package

The following packages are needed:

	1.  build-essential	(>= v11.4)
	2.  libncurses5-dev	(>= v5.7)
	3.  libpcre3-dev	(>= v7.8)
	4.  ghostscript		(>= v8.71)
	5.  libpoppler-qt4-dev	(>= v0.12)
	6.  libqt4-dev		(>= v4.6.2)


== Download ==

1. Download ChIPseeqer source code:

	http://physiology.med.cornell.edu/faculty/elemento/lab/CS_files/ChIPseeqer-2.0.tar.gz
	
2. Untar the source code archive

	tar xzf ChIPseeqer-2.0/tar.gz


== Data ==


1. Download the annotation data for a species:

	http://physiology.med.cornell.edu/faculty/elemento/lab/CS_files/
	
Untar the annotation data inside the data directory.

	cd ChIPseeqer-2.0/data/
	
	tar xzf hg18.tar.gz
	
2. Download the encode datasets:

	http://physiology.med.cornell.edu/faculty/elemento/lab/CS_files/Encode.tar.gz
	
Untar the encode datasets inside the data/NONGENIC_ANNOTATION directory.

	cd ChIPseeqer-2.0/data/NONGENIC_ANNOTATION/
	
	tar xzf Encode.tar.gz

3. To run the conservation analysis, you will need the conservation scores.

Follow these instructions to download the conservation files in your computer:

	http://icb.med.cornell.edu/wiki/index.php/Elementolab/ChIPseeqer_Cons

4. To run FIRE you will need the genome sequences and pre-packaged sequence file(s) 
for your organism.

Follow these instructions to download the required files:

	http://icb.med.cornell.edu/wiki/index.php/Elementolab/ChIPseeqer_FIRE


== How to install ==

1. Go to the source directory.

	cd ChIPseeqer-2.0/src/

2. Run the Makefile.init to retrieve information about the OS and architecture.

	make -f Makefile.init
	
3. Run make to compile basic libraries needed.

	make

4. Run make dist to create the basic tools of ChIPseeqer.

	make dist
	
5. Run make tools to create the tools used by ChIPseeqer (e.g., FIRE, PAGE).

	make tools

6. Run make gui to create the interface of ChIPseeqer.

	make gui

7. Set the envormental variables.

	Edit the .bash_profile file in your home directory and add the following at the end.
	
	export CHIPSEEQERDIR=$HOME/ChIPseeqer-2.0/dist
	export FIREDIR=$HOME/ChIPseeqer-2.0/tools/FIRE
	export PAGEDIR=$HOME/ChIPseeqer-2.0/tools/PAGE
	export HEATMAPDIR=$HOME/ChIPseeqer-2.0/tools/HEATMAP
	export PERLMODULESDIR=$HOME/ChIPseeqer-2.0/tools/PERL_MODULES
	export PERLSCRIPTSDIR=$HOME/ChIPseeqer-2.0/tools/PERL_SCRIPTS
	export MYSCANACEDIR=$HOME/ChIPseeqer-2.0/tools/MYSCANACE
	export KOHONENDIR=$HOME/ChIPseeqer-2.0/tools/KOHONEN
	
	In the example above, we assume that ChIPseeqer is installed in 
	ChIPseeqer-2.0/, inside the user's home directory. In any other case,
	make sure to edit the environment variables accordingly, so as to point 
	at the correct directories.


== How to run ChIPseeqer ==

1. To run the ChIPseeqer tools from the command-line go to the directory:

	cd ChIPseeqer-2.0/dist/
	
	Run any of the tools found at:
	
	http://icb.med.cornell.edu/wiki/index.php/Elementolab/ChIPseeqer_Tutorial

2. To run the graphical user interface of ChIPseeqer: 

	Go to the directory:

		cd ChIPseeqer-2.0/gui/

	Type:

		./ChIPseeqerGUI


== Cleanup ==
	
If needed to reinstall, go to the src directory

	cd ChIPseeqer-2.0/src/

and do the following:

	make clean
	
	make tools_clean
	
	make dist_clean
	
	make gui_clean

	make -f Makefile.init clean


== Contact ==

If you find any errors during installation, broken links, or have any comments, contact: Eugenia Giannopoulou <eug2002@med.cornell.edu>
