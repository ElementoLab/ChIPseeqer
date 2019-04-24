#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/olly\///'`
if [[ $2 != "" ]]
then ssh olivier@atlas.igh.cnrs.fr "mkdir -p ~/$dir"
fi
scp $1 olivier@atlas.igh.cnrs.fr:$dir/$1
