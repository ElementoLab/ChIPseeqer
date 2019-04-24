#!/bin/bash
file=$1
dir=`pwd | sed 's/\/home\/olly\///'`
if [[ $2 != "" ]]
then ssh website@tavazoielab.princeton.edu "mkdir -p /home/website/www/html/fastcompare/$dir"
fi
scp $1 website@tavazoielab.princeton.edu:/home/website/www/html/fastcompare/$dir/$1
