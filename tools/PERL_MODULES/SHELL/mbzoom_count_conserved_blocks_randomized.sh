# randomize
# I1 = kmer
# I2 = fasta
for r in `seq 1 1000`;
do
fasta=RANDOMIZED/DATA/$2.$r
count=RANDOMIZED/COUNTS/$1.counts.$r
/home/elemento/PROGRAMS/GENREGEXP/generate_randomized_fasta_sequences $2 > RANDOMIZED/DATA/$2.$r
mypwd=`pwd`
cat> script.pbs <<EOF
#PBS -l ncpus=1
#PBS -l mem=100Mb
#PBS -l walltime=23:00:00
#PBS -o $mypwd/ERR/$1.$r.out
#PBS -e $mypwd/ERR/$1.$r.err
cd \$PBS_O_WORKDIR
/home/elemento/PROGRAMS/GENREGEXP/genregexp_batch -refile $1 -fastafile $fasta > $count
EOF
cat script.pbs
#qsub script.pbs
rm -f $fasta
done


