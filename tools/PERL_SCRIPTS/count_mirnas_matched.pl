use strict;

my $mirnas = shift @ARGV;

my $txt = join(" ", @ARGV);

for (my $i=50; $i<=1000; $i+=50) {
    my $todo = "perl ~/PERL_MODULES/SCRIPTS/replace_7mers_by_better_scoring_kmers_single_strand.pl $i $txt > t";
    system($todo);

   my $todo = "perl ~/PERL_MODULES/SCRIPTS/fast_match_micrornas_to_kmers.pl $mirnas t 1 | columns.pl 1 | sort | uniq | wc -l";

    #my $todo = "perl ~/PERL_MODULES/SCRIPTS/match_micrornas_to_kmers.pl ~/DATA/DROSOPHILA/MIRNA/mirnas_v2.seq t 1 | columns.pl 1 | sort | uniq | wc -l";
    

    my $wc = `$todo`; 
    my ($cnt) = $wc =~ /^\ +(\d+)$/;
    
    print "$i\t$cnt\n"; 
    
}


