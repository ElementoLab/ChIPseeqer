perl mask.pl /home/olly/GENOMES/HUMAN_MOUSE/human_upstream.fasta > /home/olly/GENOMES/HUMAN_MOUSE/human_upstream_masked.fasta
perl mask.pl /home/olly/GENOMES/HUMAN_MOUSE/mouse_upstream.fasta > /home/olly/GENOMES/HUMAN_MOUSE/mouse_upstream_masked.fasta
cd /home/olly/GENOMES/HUMAN_MOUSE
/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/FASTCOMPARE/fastcompare -fasta1 human_upstream_masked.fasta -fasta2 mouse_upstream_masked.fasta -nbgenes 16810 -kmers 7mers.txt

