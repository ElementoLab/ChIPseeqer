#!/bin/bash

#
# super simple script to backup current directory to another machine
#
# this script is made much more efficient if you can access the 
# external machine without having to login ..
# for that, you have to run 
#  ssh-keygen -t rsa
#  add the newly created public key in ~/.ssh/id_rsa.pub file 
#    to ~/.ssh/authorized_keys on the other machine
#

#
# then, in any directory, just type 
#  rsynch_to_della.sh 
# to save the current directory ... (to della:/tigress/elemento/backup here)
#

machine=elemento@della.princeton.edu
todir=/tigress/elemento
dir=`pwd | sed 's/\/home\/elemento\///'`
ssh $machine "mkdir -p $todir/backup/$dir"
rsync -avz -e ssh ./ $machine:$todir/backup/$dir
