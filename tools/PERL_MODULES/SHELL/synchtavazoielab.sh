#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/elemento\///'`
if [[ $2 != "" ]]
then ssh elemento@tavazoielab.princeton.edu "mkdir -p ~/$dir"
fi
echo "scp $1 elemento@tavazoielab.princeton.edu:$dir/$1"
scp $1 elemento@tavazoielab.princeton.edu:$dir/
