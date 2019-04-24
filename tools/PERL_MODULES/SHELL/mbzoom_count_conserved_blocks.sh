rm xa*
n=`wc -l $1 | awk '{ print $1 }'`
nn=`perl -e "print int($n/12)+1"`
split -l $nn $1
mypwd=`pwd`
for f in `ls xa*`;
do
cat> script.pbs <<EOF
#PBS -l ncpus=1
#PBS -l mem=100Mb
#PBS -l walltime=10:00:00
#PBS -o $mypwd/ERR/$f.out.$c.$j
#PBS -e $mypwd/ERR/$f.err.$c.$j
cd \$PBS_O_WORKDIR
/home/elemento/PROGRAMS/GENREGEXP/genregexp_batch -refile $f -fastafile $2 > $f.counts
EOF
cat script.pbs
#qsub script.pbs
done
cat `ls x*.counts` > $1.counts

