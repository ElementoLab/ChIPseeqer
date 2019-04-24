km=$2
k7=$3
k8=$4
k9=$5
mi=$1
n7=`awk '{ if (length($1) == 7) print $0 }' $km | wc -l | sed "s/ //g"`
n8=`awk '{ if (length($1) == 8) print $0 }' $km | wc -l | sed "s/ //g"`
n9=`awk '{ if (length($1) == 9) print $0 }' $km | wc -l | sed "s/ //g"`
>RAND1
for f in `seq 1 100`;
do
>RAND
~/PERL_MODULES/SHELL/shuffle $k7 | head -$n7 >> RAND
~/PERL_MODULES/SHELL/shuffle $k8 | head -$n8 >> RAND
~/PERL_MODULES/SHELL/shuffle $k9 | head -$n9 >> RAND
perl /home/olly/PERL_MODULES/SCRIPTS/fast_match_micrornas_to_kmers.pl $mi RAND | awk '{ if ($5 < 2) print $0 }' > matches.txt
n1r=`cat matches.txt | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 0 | sort | uniq | wc -l | sed "s/ //g"`
echo $n1r
echo $n1r >> RAND1
done
n1r_avg=`perl ~/PERL_MODULES/SCRIPTS/get_average.pl RAND1`
echo " $n1 $n1r_avg"


