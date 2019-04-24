#!/bin/bash
file=$1
dir=`pwd | sed 's/\/Users\/olivier\///'`
if [[ $2 != "" ]]
then ssh elemento@della.princeton.edu "mkdir -p ~/$dir"
fi
echo "scp $1 elemento@della.princeton.edu:$dir/$1"
scp $1 elemento@della.princeton.edu:$dir/$1
