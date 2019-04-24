echo "Number of sequences :"
wc -l $1
getseqsfromblastdb.pl $2 < $1 > tmp
/home/olly/CELEGANS/ALIGNACE/AlignACE -i tmp -numcols 10 -oversample 2 | tee AlignACE.out

