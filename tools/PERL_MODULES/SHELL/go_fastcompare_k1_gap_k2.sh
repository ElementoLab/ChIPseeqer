for k in `seq 6 8`;
do
file1=`echo -n $k && echo -n "mers.txt"`
for k1 in `seq 2 4`;
do
for gap in `seq 2 4`;
do
/home/olly/PROGRAMS/FASTCOMPARE/fastcompare_fast_old_single_strand -kmers /home/olly/DATA/KMERS/ALL/$file1 -fasta1 ../DATA/ELE_D_300.seq -fasta2 ../DATA/BRI_D_300.seq -nbgenes 11292 -k $k -k1 $k1 -gap $gap -singlestrand T >  $k.$k1.$gap.txt
done
done
done
