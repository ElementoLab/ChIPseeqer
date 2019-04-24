#!/bin/bash
# 6, 7, 8, 9 and 10 mers
data1=$1
data2=$2
resdir=$4
nb=$3
lists=/home/olly/DATA/KMERS/FILTERED
#options="-singlestrand Y"
#nb=4358
for i in `seq 7 10`;
do
echo -n $i && echo "mers .."
file1=`echo -n $i && echo -n "mers.txt"`
file2=`echo -n $i && echo -n "mers_all.txt"`
echo "./fastcompare_fast -k $i -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/$file1 $options > $resdir/$file2"
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand -k $i -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/$file1 $options > $resdir/$file2
done
# 6 mers, gaps 1->10
for i in `seq 2 20`;
do
echo -n "6mers, $i bp gap .."
file2=`echo -n "6mers_gap" && echo -n $i && echo -n "all.txt"`
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand -gap $i -k 6 -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/6mers.txt $options > $resdir/$file2
done
# 8 mers, gaps 1->10
for i in `seq 2 20`;
do
echo -n "8mers, $i bp gap .."
file2=`echo -n "8mers_gap" && echo -n $i && echo -n "all.txt"`
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand -gap $i -k 8 -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/8mers.txt $options > $resdir/$file2
done
for i in `seq 2 20`;
do
echo -n "10mers, $i bp gap .."
file2=`echo -n "10mers_gap" && echo -n $i && echo -n "all.txt"`
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand -gap $i -k 10 -fasta1 $data1 -fasta2 $data2 -nbgenes $nb -kmers $lists/10mers.txt $options > $resdir/$file2
done

