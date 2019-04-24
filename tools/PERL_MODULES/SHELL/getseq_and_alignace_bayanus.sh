echo "Number of sequences :"
wc -l $1
sh /home/olly/PERL_MODULES/bayanus_upstream_regions.sh $1 > tmp
/home/olly/CELEGANS/ALIGNACE/AlignACE -i tmp | tee AlignACE.out

