# input : fasta file, homologie file

# generate FIRE profile input
shuffle_fasta_file_markov.pl $1 > $1.shu
fasta_sequence_names.pl $1.shu | awk '{ print $1 "\t" 0 }' > $1.exp.c0
fasta_sequence_names.pl $1     | awk '{ print $1 "\t" 1 }' > $1.exp.c1
cat <<EOF > h.exp
i	i
EOF
cat h.exp $1.exp.c0 $1.exp.c1 > $1.exp 
echo generated $1.exp

# generate homology
fasta_sequence_names.pl $1.shu > $1.shu.homologies
cat $2 $1.shu.homologies > $1.all.homologies
echo generated $1.all.homologies


# generate fasta file
cat $1 $1.shu > $1.all
echo generated $1.all

# menage
\rm h.exp $1.exp.c0 $1.exp.c1
\rm $1.shu
\rm $1.shu.homologies