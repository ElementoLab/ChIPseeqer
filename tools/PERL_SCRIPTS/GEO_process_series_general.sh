# usage: sh ../GEO_process_series.sh GSE2836_series_matrix.txt GPL2569.txt 0 8 genelist delim

scriptdir=/home/elemento/PERL_MODULES/SCRIPTS

#
# step 1 : clean series matrix
#
echo "1. cleaning matrix"
perl -pi -e 's/\r//g' $1
perl -pi -e 's/\r//g' $2
perl $scriptdir/GEO_series_to_matrix.pl $1 > $1.clean


#
# step 2 : clean platform, get mapping
#
echo "2. cleaning platform, getting raw mapping"
perl ~/PERL_MODULES/SCRIPTS/columns.pl $3 $4 < $2  > $2.clean

#
# step 3 : cleanup mapping
#
echo "3 cleanup mapping"
perl ~/PERL_MODULES/SCRIPTS/ids_consolidate_mapping_using_gene_list.pl $2.clean $5 > $2.clean.simple
cp $2.clean.simple $2.clean

# perl ~/PERL_MODULES/SCRIPTS/ids_consolidate_tabdelim_mapping_using_gene_list.pl $2.clean /home/elemento/PROGRAMS/FACE/TEST/ANNOTATIONS/human_names_from_go_index.txt > $2.clean.simple
# perl ~/PERL_MODULES/SCRIPTS/ids_consolidate_tabdelim_mapping_using_gene_list.pl $2.clean $tmpdir/mouse_go_index.txt.genenames > $2.clean.simple
# perl ~/PERL_MODULES/SCRIPTS/ids_consolidate_commadelim_mapping_using_gene_list.pl $2.clean /home/elemento/PROGRAMS/FACE/TEST/ANNOTATIONS/human_names_from_go_index.txt > $2.clean.simple


echo "4. translate ids"
perl  $scriptdir/translate_column_using_table_column.pl --table=$1.clean --col=0 --dict=$2.clean --k=0 --v=1 > $1.ids


echo "5. average rows"
perl ~/PERL_MODULES/SCRIPTS/expression_average_rows_with_same_id.pl $1.ids > $1.rowavg
wc -l $1.rowavg
#perl $tmpdir/expression_average_all_columns.pl $1.rowavg > $1.final
#exit
#\rm $1.ids
#\rm $2.clean
#\rm $1.clean
\rm *.bak
exit

echo "5. normalize"
#perl ~/PERL_MODULES/SCRIPTS/expression_normalize_rows.pl < $1.rowavg > $1.final
a=`wc -l $1.final | awk '{ print $1 }'`
n=`perl -e "print int(0.5+sqrt($a))"`
echo "6. cluster into $n groups"
~/PROGRAMS/MIMOTIFS/TEMPORARY/CLUSTER/cluster-1.35/src/cluster.exe -k $n -r 10 -f $1.final
