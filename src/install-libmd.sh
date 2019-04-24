#!/bin/sh
if [ ! -e libmd/libmd.a ]
then
	echo "installing libmd"
	(cd libmd; perl Makefile.PL; make; cd ..)
fi 

