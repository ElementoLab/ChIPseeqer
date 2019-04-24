mysql WORMS -e "select * from GO_FULL_ANNOTATION" | awk '{ print $2 }' | sort | uniq > tmp1
fasta_sequence_names.pl ../DATA/ELE_D_300.seq > tmp2
getUnion.pl tmp1 tmp2 | wc -l
rm tmp1
rm tmp2
