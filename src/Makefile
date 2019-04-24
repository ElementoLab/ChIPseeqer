# typical flags and options
CFLAGS	= -Wall -O1
LDFLAGS	=
CC	= cc
C++	= c++
DIST	= ../dist
TOOLS	= ../tools
GUI	= ../gui
SRC	= ../src
INCC	= 

# DO NOT EDIT BELOW THIS LINE

CFLAGS	+= -I lib/third_party/zlib/
LDFLAGS	+= -lm -Llib/third_party/zlib/	\
	lib/third_party/samtools/libbam.a lib/third_party/zlib/libz.a

# architecture and OS detection 
OS 	= $(shell ./os)
ARCH	= $(shell ./arch)

# debug 
ifeq (1, $(DEBUG))
CFLAGS 	+= -g -ggdb
INCC	+= -g -ggdb
endif

# OS
ifeq ($(OS), darwin)
# arch
ifeq ($(ARCH), i386)
CFLAGS	+= -arch i386
INCC	+= -arch i386
endif
ifeq ($(ARCH), x86_64)
CFLAGS	+= -arch x86_64
INCC	+= -arch x86_64
endif
ifeq (1, $(shell [ -d /opt/local64 ] && echo 1 ))
CFLAGS	+= -I/opt/local64/include 
LDFLAGS	+= -L/opt/local64/lib
INCC	+= -I/opt/local64/include
INCL	= -L/opt/local64/lib
endif
ifeq (1, $(shell [ -d /opt/local ] && echo 1 ))
CFLAGS	+= -I/opt/local/include 
LDFLAGS	+= -L/opt/local/lib
INCC	+= -I/opt/local/include
INCL	= -L/opt/local/lib
endif
endif

PROGCHIP		= ChIPseeqer.bin
PROGCONS		= ChIPseeqer2Cons.bin
PROGCLONAL		= ChIPseeqerGetNumClonalReads
PROGFIND		= ChIPseeqerGetReadCountInPeakRegions
PROGDENS		= ChIPseeqerGetReadDensityProfiles.bin
PROGGENEDATA    = ChIPseeqerMakeGeneEnrichmentData
PROGCNTDIST		= ChIPseeqerReadCountDistribution
PROGCMP			= CompareIntervals
PROGDIST		= FindClosestGene
PROGGENFRAC     = FindGenomeFractions
PROGGENREGEXP	= genregexp
PROGSUB			= GetSubstrings
PROGSUBIDX		= GetSubstringsIdx
PROGRANDREG		= GetRandomRegions
PROGRANDREGCPG	= GetRandomRegionsCpG
PROGRANDREGFILE = GetRandomRegionsFromFile
PROGFASTAIDX    = IndexFasta
PROGDENSTRACK	= MakeGenomicReadDensityTrack
PROGMAP			= MakeMappabilityStats
PROGMAPMAP		= MakeMappabilityChrMaps
PROGHYPERGEOM	= myhypergeom
PROGEXTRAND		= ExtractRandomLines
PROGSPLSAM		= split_samfile
PROGSPLBAM		= split_bamfile

.PHONY: all dist clean dist_clean tools tools_clean gui gui_clean

all : $(PROGFIND) $(PROGCHIP) $(PROGCMP) $(PROGDIST) $(PROGSUB)		\
	$(PROGCONS) $(PROGDENS) $(PROGMAP) $(PROGSPLSAM)		\
	$(PROGGENEDATA) $(PROGDENSTRACK) $(PROGMAPMAP)	\
	$(PROGGENREGEXP) $(PROGSUBIDX) $(PROGCLONAL) $(PROGSPLBAM)      \
	$(PROGFASTAIDX) $(PROGGENFRAC) $(PROGRANDREG) $(PROGEXTRAND)	\
	$(PROGHYPERGEOM) $(PROGCNTDIST) $(PROGRANDREGCPG) $(PROGRANDREGFILE)

$(PROGGENEDATA): $(PROGGENEDATA).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGGENEDATA) $(PROGGENEDATA).c		\
		statistics.o information.o dataio.o hashtable.o set.o	\
		sequences.o  readio.o $(LDFLAGS)

$(PROGGENFRAC): $(PROGGENFRAC).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGGENFRAC) $(PROGGENFRAC).c		\
		statistics.o information.o dataio.o hashtable.o set.o	\
		sequences.o  readio.o $(LDFLAGS)

$(PROGCHIP): $(PROGCHIP).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o libmd
	@sh install-libmd.sh
	$(CC) $(CFLAGS) -o $(PROGCHIP) $(PROGCHIP).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o  libmd/libmd.a $(LDFLAGS)

$(PROGPROF): $(PROGPROF).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGPROF) $(PROGPROF).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGDENS): $(PROGDENS).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o postscript.o readio.o
	$(CC) $(CFLAGS) -o $(PROGDENS) $(PROGDENS).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		postscript.o readio.o $(LDFLAGS)

$(PROGFIND): $(PROGFIND).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGFIND) $(PROGFIND).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGCMP): $(PROGCMP).cc statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o interval_tree.o readio.o
	$(C++) $(CFLAGS) -o $(PROGCMP) $(PROGCMP).cc statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		interval_tree.o readio.o $(LDFLAGS)

$(PROGDIST): $(PROGDIST).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGDIST) $(PROGDIST).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGCONS): $(PROGCONS).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o postscript.o
	$(CC) $(CFLAGS) -o $(PROGCONS) $(PROGCONS).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		postscript.o readio.o -lz $(LDFLAGS)

$(PROGSUB): $(PROGSUB).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGSUB) $(PROGSUB).c statistics.o		\
		information.o dataio.o hashtable.o readio.o set.o	\
		sequences.o $(LDFLAGS)

$(PROGSUBIDX): $(PROGSUBIDX).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGSUBIDX) $(PROGSUBIDX).c statistics.o	\
		information.o dataio.o hashtable.o readio.o set.o	\
		sequences.o $(LDFLAGS)
		 
$(PROGRANDREG): $(PROGRANDREG).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGRANDREG) $(PROGRANDREG).c statistics.o	\
		information.o dataio.o hashtable.o readio.o set.o	\
		sequences.o $(LDFLAGS)
		
$(PROGRANDREGCPG): $(PROGRANDREGCPG).cc statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o interval_tree.o readio.o
	$(C++) $(CFLAGS) -o $(PROGRANDREGCPG) $(PROGRANDREGCPG).cc statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		interval_tree.o readio.o $(LDFLAGS)
		
$(PROGRANDREGFILE): $(PROGRANDREGFILE).cc statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o interval_tree.o readio.o
	$(C++) $(CFLAGS) -o $(PROGRANDREGFILE) $(PROGRANDREGFILE).cc statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		interval_tree.o readio.o $(LDFLAGS)

$(PROGMAP): $(PROGMAP).c statistics.o information.o dataio.o		\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGMAP) $(PROGMAP).c statistics.o		\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGMAPMAP): $(PROGMAPMAP).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGMAPMAP) $(PROGMAPMAP).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGSPLSAM): $(PROGSPLSAM).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGSPLSAM) $(PROGSPLSAM).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGSPLBAM): $(PROGSPLBAM).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGSPLBAM) $(PROGSPLBAM).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGDENSTRACK): $(PROGDENSTRACK).c statistics.o information.o		\
	dataio.o hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGDENSTRACK) $(PROGDENSTRACK).c		\
		statistics.o information.o dataio.o hashtable.o set.o	\
		sequences.o readio.o $(LDFLAGS)

$(PROGCNTDIST): $(PROGCNTDIST).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGCNTDIST) $(PROGCNTDIST).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGCLONAL): $(PROGCLONAL).c statistics.o information.o dataio.o	\
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS) -o $(PROGCLONAL) $(PROGCLONAL).c statistics.o	\
		information.o dataio.o hashtable.o set.o sequences.o	\
		readio.o $(LDFLAGS)

$(PROGFASTAIDX): $(PROGFASTAIDX).c statistics.o information.o dataio.o  \
	hashtable.o set.o sequences.o readio.o
	$(CC) $(CFLAGS)  -o $(PROGFASTAIDX) $(PROGFASTAIDX).c 		\
		statistics.o information.o dataio.o hashtable.o set.o   \
		sequences.o readio.o $(LDFLAGS)

$(PROGGENREGEXP): $(PROGGENREGEXP).c
	$(CC) $(CFLAGS) -lpcre -o $(PROGGENREGEXP) $(PROGGENREGEXP).c $(LDFLAGS)
	
$(PROGEXTRAND): $(PROGEXTRAND).c
	$(CC) $(CFLAGS) -o $(PROGEXTRAND) $(PROGEXTRAND).c $(LDFLAGS)

$(PROGHYPERGEOM): $(PROGHYPERGEOM).c
	$(CC) $(CFLAGS) -o $(PROGHYPERGEOM) $(PROGHYPERGEOM).c $(LDFLAGS)

interval_tree.o: interval_tree.c interval_tree.h templatestack.h
	$(C++) $(CFLAGS) -c interval_tree.c 

mi_library.o: mi_library.c mi_library.h
	$(CC) $(CFLAGS) -c mi_library.c 

statistics.o: statistics.c statistics.h
	$(CC) $(CFLAGS) -c statistics.c 

sequences.o: sequences.c sequences.h
	@sh install.sh
	$(CC) $(CFLAGS) -c sequences.c 

dataio.o: dataio.c dataio.h
	$(CC) $(CFLAGS) -c dataio.c 

readio.o: readio.c readio.h
	$(CC) $(CFLAGS) -c readio.c 

information.o: information.c information.h
	$(CC) $(CFLAGS) -c information.c 

regexp.o: regexp.c regexp.h
	$(CC) $(CFLAGS) -c regexp.c 

prefix.o: prefix.c prefix.h
	$(CC) $(CFLAGS) -c prefix.c 

hashtable.o: hashtable.c hashtable.h
	$(CC) $(CFLAGS) -c hashtable.c

set.o: set.c set.h
	$(CC) $(CFLAGS) -c set.c

postscript.o: postscript.c postscript.h
	$(CC) $(CFLAGS) -c postscript.c  

dist:
	mv $(PROGCHIP)		$(DIST)
	mv $(PROGCONS)		$(DIST)
	mv $(PROGCLONAL)	$(DIST)
	mv $(PROGFIND)		$(DIST)
	mv $(PROGDENS)		$(DIST)
	mv $(PROGGENEDATA)	$(DIST)
	mv $(PROGCNTDIST)	$(DIST)
	mv $(PROGCMP)		$(DIST)
	mv $(PROGDIST)		$(DIST)
	mv $(PROGGENFRAC)	$(DIST)
	mv $(PROGGENREGEXP)	$(DIST)
	mv $(PROGSUB)		$(DIST)
	mv $(PROGSUBIDX)	$(DIST)
	mv $(PROGRANDREG)	$(DIST)
	mv $(PROGRANDREGCPG)	$(DIST)
	mv $(PROGRANDREGFILE)	$(DIST)
	mv $(PROGFASTAIDX)	$(DIST)
	mv $(PROGDENSTRACK)	$(DIST)
	mv $(PROGMAP)		$(DIST)
	mv $(PROGMAPMAP)	$(DIST)
	mv $(PROGHYPERGEOM) $(DIST)
	mv $(PROGEXTRAND)	$(DIST)
	mv $(PROGSPLSAM)	$(DIST)
	mv $(PROGSPLBAM)	$(DIST)

clean:
	@sh uninstall.sh
	@sh uninstall-libmd.sh
	rm -rf	*.dSYM *.o $(PROGFIND) $(PROGCHIP) $(PROGCMP) $(PROGDIST)	\
		$(PROGSUB) $(PROGCONS) $(PROGPROF) $(PROGGENREGEXP)		\
		$(PROGCNTDIST) $(PROGCLONAL) $(PROGDENS)			\
		$(PROGGENEDATA) $(PROGSUBIDX) $(PROGDENSTRACK)			\
		$(PROGMAPMAP) $(PROGMAP) $(PROGSPLSAM)	$(PROGSPLBAM)		\
		$(PROGRANDREG) $(PROGEXTRAND) $(PROGGENFRAC) $(PROGFASTAIDX)	\
		$(PROGRANDREGCPG) $(PROGRANDREGFILE) *.pro *.xcodeproj

dist_clean:
	rm -f	$(DIST)/$(PROGFIND) $(DIST)/$(PROGCHIP)		\
		$(DIST)/$(PROGCMP) $(DIST)/$(PROGDIST)		\
		$(DIST)/$(PROGSUB) $(DIST)/$(PROGCONS)		\
		$(DIST)/$(PROGGENREGEXP)			\
		$(DIST)/$(PROGCNTDIST) $(DIST)/$(PROGCLONAL)	\
		$(DIST)/$(PROGDENS) $(DIST)/$(PROGGENEDATA)	\
		$(DIST)/$(PROGSUBIDX) $(DIST)/$(PROGDENSTRACK)	\
		$(DIST)/$(PROGMAPMAP) $(DIST)/$(PROGMAP)	\
		$(DIST)/$(PROGSPLSAM) $(DIST)/$(PROGSPLBAM)	\
		$(DIST)/$(PROGFASTAIDX) $(DIST)/$(PROGGENFRAC) 	\
		$(DIST)/$(PROGRANDREG) $(DIST)/$(PROGHYPERGEOM) \
		$(DIST)/$(PROGEXTRAND) $(DIST)/$(PROGRANDREGCPG) \
		$(DIST)/$(PROGRANDREGFILE)

tools:
	cd $(TOOLS)/FIRE && (CFLAGS="$(INCC)" LDFLAGS="$(INCL)" ./configure) && cd ..
	make -C $(TOOLS)/FIRE
	make -C $(TOOLS)/PAGE CFLAGS="$(INCC)"
	make -C $(TOOLS)/MYSCANACE CFLAGS="$(INCC)"
	make -C $(TOOLS)/KOHONEN CFLAGS="$(INCC)"

tools_clean:
	make -C $(TOOLS)/FIRE clean
	make -C $(TOOLS)/PAGE clean
	make -C $(TOOLS)/MYSCANACE clean
	make -C $(TOOLS)/KOHONEN clean

gui:
ifeq ($(OS), darwin)
	cp $(GUI)/commonPaths.h.tmpl $(GUI)/commonPaths.h
endif
ifeq ($(OS), linux)
	cat $(GUI)/commonPaths.h.tmpl | sed -e 's/\"\.\.\/\.\.\/\.\.\//"/g' > $(GUI)/commonPaths.h
endif
	echo "cd $(GUI)" > build.gui
ifeq ($(ARCH), i386)
	echo "cat ChIPseeqerGUI.pro | sed -e 's/^CONFIG\(.\)*/CONFIG += x86/g' > ChIPseeqerGUI.pro.tmp" >> build.gui
endif
ifeq ($(ARCH), x86_64)
	echo "cat ChIPseeqerGUI.pro | sed -e 's/^CONFIG\(.\)*/CONFIG += x86_64/g' > ChIPseeqerGUI.pro.tmp" >> build.gui
endif
	echo "mv ChIPseeqerGUI.pro.tmp ChIPseeqerGUI.pro" >> build.gui
	echo "qmake" >> build.gui
	echo "cd $(SRC)" >> build.gui
	chmod +x build.gui
	./build.gui
	make -C $(GUI)

gui_clean:
	make -C $(GUI) clean
ifeq ($(OS), linux)
	rm -rf $(GUI)/ChIPseeqerGUI
endif
	rm -rf $(GUI)/commonPaths.h
	rm -rf $(GUI)/Makefile
	rm -rf build.gui
