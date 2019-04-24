#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/elemento\/PROGRAMS\/MIMOTIFS\//www\/html\/FIRE\//'`
if [[ $2 != "" ]]
then ssh website@tavazoielab.princeton.edu "mkdir -p ~/$dir"
fi
echo "scp $1 website@tavazoielab.princeton.edu:$dir/$1"
scp $1 website@tavazoielab.princeton.edu:$dir/
