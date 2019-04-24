#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/olly\///'`
if [[ $2 != "" ]]
then ssh elemento@straylite "mkdir -p ~/$dir"
fi
scp $1 elemento@straylite:$dir/$1
