#!/bin/sh
# uninstall third-party libraries
#	1. samtools
#	2. zlib

echo "uninstalling samtools"
	(cd lib/third_party/samtools && make clean)

echo "uninstalling zLib"
	(cd lib/third_party/zlib-1.2.5 && make clean && rm -Rf build_dir && cd ..; rm -rf zlib)
