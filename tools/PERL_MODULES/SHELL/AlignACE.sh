#!/bin/bash
for c in `seq 6 12`;
do
for j in `seq 1 5`;
do
for f in `ls *.fsa`; 
do
echo "/home/olly/PERL_MODULES/PROGRAMS/ACE/AlignACE -numcols $c -i $f -oversample 2 > OUT/$f.out.$c.$j"
/home/olly/PERL_MODULES/PROGRAMS/ACE/AlignACE -numcols $c -i $f -oversample 2 > OUT/$f.out.$c.$j
done
done
done
