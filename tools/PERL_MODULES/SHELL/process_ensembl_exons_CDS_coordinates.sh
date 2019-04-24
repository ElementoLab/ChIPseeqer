vim $1 -s /home/olly/PERL_MODULES/SCRIPTS/vi_dell_first_line.txt
columns.pl 2 0 5 6 1 3 4 < $1 > $1-rearranged
perl ~/PERL_MODULES/SCRIPTS/extract_gene_boundaries_from_exon_and_CDS_boundaries.pl $1-rearranged > $1-rearranged-with_3UTRs
perl ~/PERL_MODULES/SCRIPTS/calculate_3UTR_lengths_from_gene_boundaries.pl $1-rearranged-with_3UTRs | awk '{ if ($2 > 2) print $0 }'  > $1-rearranged-with_3UTRs-lengths
