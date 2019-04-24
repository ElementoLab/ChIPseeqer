kmers=$1
mirnas=$2
k7=$3
k8=$4
k9=$5
perl /home/olly/PERL_MODULES/SCRIPTS/match_micrornas_to_kmers.pl $2 $1 | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 3 | tee pos.txt
n7=`awk '{ if (length($1) == 7) print $0 }' $1 | wc -l | sed "s/ //g"`
n8=`awk '{ if (length($1) == 8) print $0 }' $1 | wc -l | sed "s/ //g"`
n9=`awk '{ if (length($1) == 9) print $0 }' $1 | wc -l | sed "s/ //g"`
for f in `seq 1 100`;
do
>RAND
perl ~/PERL_MODULES/SCRIPTS/shuffle_rows.pl $k7 | head -$n7 >> RAND
perl ~/PERL_MODULES/SCRIPTS/shuffle_rows.pl $k8 | head -$n8 >> RAND
perl ~/PERL_MODULES/SCRIPTS/shuffle_rows.pl $k9 | head -$n9 >> RAND
perl /home/olly/PERL_MODULES/SCRIPTS/match_micrornas_to_kmers.pl $2 RAND | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 3 | tee pos$f.txt
done

