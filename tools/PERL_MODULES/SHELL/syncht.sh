#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/elemento\///'`
if [[ $2 != "" ]]
then ssh elemento@tcluster.princeton.edu "mkdir -p ~/$dir"
fi
echo "scp $1 elemento@tcluster.princeton.edu:$dir/"
scp $1 elemento@tcluster.princeton.edu:$dir/
