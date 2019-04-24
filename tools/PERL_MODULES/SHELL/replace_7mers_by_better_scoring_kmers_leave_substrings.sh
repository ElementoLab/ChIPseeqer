#!/bin/bash
perl /home/olly/PROGRAMS/CLUSTERING_KMERS/replace_7mers_by_better_scoring_kmers.pl  $1 $2 > /tmp/tmp-clu.txt
perl /home/olly/PROGRAMS/CLUSTERING_KMERS/replace_7mers_by_better_scoring_kmers.pl /tmp/tmp-clu.txt $3 | sort | uniq | perl ~/PERL_MODULES/SCRIPTS/sort_column_inv.pl 4 | perl ~/PERL_MODULES/SCRIPTS/put_inf_first.pl
rm /tmp/tmp-clu.txt

