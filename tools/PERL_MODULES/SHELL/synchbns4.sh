#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/elemento\///'`
if [[ $2 != "" ]]
then ssh olly@bns4.princeton.edu "mkdir -p ~/$dir"
fi
echo "scp $1 olly@bns4.princeton.edu:$dir/$1"
scp $1 olly@bns4.princeton.edu:$dir/$1
