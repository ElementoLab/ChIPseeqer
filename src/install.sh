#!/bin/sh
# install third-party libraries
#	1. samtools
#	2. zlib

if [ ! -f lib/third_party/zlib/libz.a ]
then
	ARCH=`./arch`
	if [ $ARCH = "x86_64" ]; then
		OPT="--64"
	else
		OPT="--32"
	fi
	echo "installing zLib"
	(cd lib/third_party/zlib-1.2.5 && mkdir build_dir && ./configure $OPT 	\
	--prefix=build_dir --libdir=build_dir		\
	--includedir=build_dir && make && make install	\
	&& cd .. && rm -rf zlib && ln -fs zlib-1.2.5/build_dir zlib)
fi 

if [ ! -f lib/third_party/samtools/libbam.a ]
then
	echo "installing samtools"
	ARCH=`./arch`
	if [ $ARCH = "x86_64" ]; then
		OPT="-m64"
	else
		OPT="-m32"
	fi
	(cd lib/third_party/samtools && make CFLAGS="$OPT")
fi 
