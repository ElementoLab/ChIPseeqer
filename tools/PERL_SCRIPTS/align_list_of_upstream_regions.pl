
use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;

my $a_ref_cg = Sets::readSet($ARGV[0]);

foreach my $r (@$a_ref_cg) {
    print "--> $r\n";
    system("perl ~/PERL_MODULES/SCRIPTS/align_upstream_or_downstream_regions.pl $r U 2000 DATA/ortholog_table.txt");
    
    system("perl ~/PERL_MODULES/SCRIPTS/extract_fasta_from_ali.pl < $r.ali > $ARGV[1]/$r.ali.fa");
    system("perl ~/PERL_MODULES/SCRIPTS/dialign2clustalw.pl $ARGV[1]/$r.ali.fa");
}
