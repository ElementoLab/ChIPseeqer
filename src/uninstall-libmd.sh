#!/bin/sh
if [ -e libmd/libmd.a ]
then
	echo "uninstalling libmd"
	(cd libmd; make clean; cd ..)
fi 

