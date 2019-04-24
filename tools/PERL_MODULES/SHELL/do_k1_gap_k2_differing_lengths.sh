#!/bin/bash
for k in `seq 6 8`;
do
file1=`echo -n $k && echo -n "mers.txt"`
for k1 in `seq 2 4`;
do
for gap in `seq 2 4`;
do
echo $k.$k1.$gap.txt
echo "/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand_differing_lengths -kmers /home/olly/DATA/KMERS/ALL/$file1 -fasta1 $1 -fasta2 $2 -nbgenes $3 -k $k -k1 $k1 -gap $gap -singlestrand T -norm $4"
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand_differing_lengths -kmers /home/olly/DATA/KMERS/ALL/$file1 -fasta1 $1 -fasta2 $2 -nbgenes $3 -k $k -k1 $k1 -gap $gap -singlestrand T -norm $4 > $k.$k1.$gap.txt
head -10 $k.$k1.$gap.txt
done
done
done
