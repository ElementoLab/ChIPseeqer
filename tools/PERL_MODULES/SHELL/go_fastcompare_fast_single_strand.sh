#!/bin/bash
# 6, 7, 8, 9 and 10 mers
data1=$1
data2=$2
resdir=$4
nb=$3
norm=$5
lists=/home/olly/DATA/KMERS/ALL
options="-singlestrand T"
#nb=4358
for i in `seq 7 7`;
do
echo -n $i && echo "mers .."
file1=`echo -n $i && echo -n "mers.txt"`
file2=`echo -n $i && echo -n "mers_all.txt"`
echo "/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand_differing_lengths -k $i -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/$file1 $options -norm $norm > $resdir/$file2"
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand_differing_lengths -k $i -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/$file1 $options -norm $norm > $resdir/$file2
done

