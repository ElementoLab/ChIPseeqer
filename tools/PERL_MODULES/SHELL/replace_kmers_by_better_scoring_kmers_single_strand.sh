#!/bin/bash
perl /home/olly/PERL_MODULES/SCRIPTS/replace_7mers_by_better_scoring_kmers.pl  $1 $2 > /tmp/tmp-clu.txt
perl /home/olly/PROGRAMS/CLUSTERING_KMERS/replace_7mers_by_better_scoring_kmers.pl /tmp/tmp-clu.txt $3 | sort | uniq | perl ~/PERL_MODULES/SCRIPTS/sort_column_inv.pl 4 | perl ~/PERL_MODULES/SCRIPTS/put_inf_first.pl > tmp1.txt && perl ~/PERL_MODULES/SCRIPTS/remove_substrings_from_kmer_list.pl tmp1.txt
rm /tmp/tmp-clu.txt
rm  tmp1.txt

