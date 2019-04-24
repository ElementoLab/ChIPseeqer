n=`wc -l $2 | awk '{ print $1 }'`
nn=`perl -e "print int($n/2)+1"`
split -l $nn $2
t1a="perl /home/olly/PROGRAMS/FASTCOMPARE/recompare_allwindows_allkmers.pl --kmerfile=xaa --species=$1"
t1b="perl /home/olly/PROGRAMS/FASTCOMPARE/recompare_allwindows_allkmers.pl --kmerfile=xab --species=$1"
t2="perl /home/olly/PERL_MODULES/SCRIPTS/evaluate_kmers.pl --species=$1 --kmerfile=$2 --dist=1 --cond=1 --tfac=1"
t3="perl /home/olly/PERL_MODULES/SCRIPTS/evaluate_kmers.pl --species=$1 --kmerfile=$2 --orie=1 --func=1 --chip=1 --mots=1 --bestwins=$2.bestwindows"
t4="perl /home/olly/PERL_MODULES/SCRIPTS/locate_kmers.pl --utr5file=$3 --kmerfile=$2 --bestwins=$2.bestwindows --orfsfile=$4"
echo $t1a
$t1a > xaa.bestwindows &
echo $t1b
$t1b > xab.bestwindows &
wait
cat xaa.bestwindows xab.bestwindows > $2.bestwindows
echo $t2
$t2 > $2.evaluation.1
echo $t3
$t3 > $2.evaluation.2
echo $t4
$t4 > $2.evaluation.4
