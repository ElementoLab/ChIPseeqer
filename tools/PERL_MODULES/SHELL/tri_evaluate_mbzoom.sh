rm xa*
n=`wc -l $2 | awk '{ print $1 }'`
nn=`perl -e "print int($n/2)+1"`
split -l $nn $2
mypwd=`pwd`
for f in `ls xa`;
do
cat> script.pbs <<EOF
#PBS -l ncpus=1
#PBS -l mem=100Mb
#PBS -l walltime=10:00:00
#PBS -o $mypwd/ERR/$f.out.$c.$j
#PBS -e $mypwd/ERR/$f.err.$c.$j
cd \$PBS_O_WORKDIR
perl /home/elemento/PROGRAMS/FASTCOMPARE/recompare_allwindows_allkmers.pl --kmerfile=$f --species=$1 > $mypwd/$f.bestwindows; 
EOF
cat script.pbs
#qsub script.pbs
done
cat `ls x*.bestwindows` > $2.bestwindows 

