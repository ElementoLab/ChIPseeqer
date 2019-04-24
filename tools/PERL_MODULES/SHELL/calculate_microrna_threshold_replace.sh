for t in `seq 50 50 1050`;
do
k7=$2
k8=$3
k9=$4
mi=$1 
echo -n $t
perl ~/PERL_MODULES/SCRIPTS/replace_7mers_by_better_scoring_kmers_single_strand.pl $t $k7 $k8 $k9 > tmp
perl /home/olly/PERL_MODULES/SCRIPTS/fast_match_micrornas_to_kmers.pl $mi tmp > matches.txt
n1=`cat matches.txt | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 1 | sort | uniq | wc -l | sed "s/ //g"`
n2=`cat matches.txt | awk '{ if ($5 < 2) print $0 }' | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 1 | sort | uniq | wc -l | sed "s/ //g"`
#echo $n1 $n2
n7=`awk '{ if (length($1) == 7) print $0 }' tmp | wc -l | sed "s/ //g"`
n8=`awk '{ if (length($1) == 8) print $0 }' tmp | wc -l | sed "s/ //g"`
n9=`awk '{ if (length($1) == 9) print $0 }' tmp | wc -l | sed "s/ //g"`
>RAND1
>RAND2
for f in `seq 1 100`;
do
>RAND
~/PERL_MODULES/SHELL/shuffle $k7 | head -$n7 >> RAND
~/PERL_MODULES/SHELL/shuffle $k8 | head -$n8 >> RAND
~/PERL_MODULES/SHELL/shuffle $k9 | head -$n9 >> RAND
perl /home/olly/PERL_MODULES/SCRIPTS/fast_match_micrornas_to_kmers.pl $mi RAND > matches.txt
n1r=`cat matches.txt | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 1 | sort | uniq | wc -l | sed "s/ //g"`
n2r=`cat matches.txt | awk '{ if ($5 < 2) print $0 }' | perl /home/olly/PERL_MODULES/SCRIPTS/columns.pl 1 | sort | uniq | wc -l | sed "s/ //g"`
echo $n1r >> RAND1
echo $n2r >> RAND2
done
n1r_avg=`perl ~/PERL_MODULES/SCRIPTS/get_average.pl RAND1`
n2r_avg=`perl ~/PERL_MODULES/SCRIPTS/get_average.pl RAND2`
echo " $n1 $n2 $n1r_avg $n2r_avg"
done

