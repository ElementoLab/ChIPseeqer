#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/olly\///'`
scp $1 elemento@mbzoom.princeton.edu:$dir/$1
