
use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;

my $a_ref_cg = Sets::readSet($ARGV[0]);

foreach my $r (@$a_ref_cg) {
    print "--> $r\n";
    #system("perl ~/PERL_MODULES/SCRIPTS/align_3utr.pl $r");
    
    #system("perl ~/PERL_MODULES/SCRIPTS/extract_fasta_from_ali.pl < UTR.ali > UTR.ali.fa");
    
    system("perl ~/PERL_MODULES/SCRIPTS/count_conserved_kmers_in_fasta_sequences.pl UTR.ali.fa");

    #system("mv UTR.ali.fa 3UTR_ALIGNMENTS/$r.ali.fa");
}
