#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/elemento\///'`
if [[ $2 != "" ]]
then ssh elemento@gen-comp2.princeton.edu "mkdir -p /Genomics/fafner/grid/users/elemento/$dir"
fi
echo "scp $1 elemento@gen-comp2.princeton.edu:/Genomics/fafner/grid/users/elemento/$dir/$1"
scp $1 elemento@gen-comp2.princeton.edu:/Genomics/fafner/grid/users/elemento/$dir/$1
