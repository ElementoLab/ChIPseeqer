use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {

    my $todo = "/home/olly/PROGRAMS/GENREGEXP/genregexp -re $r->[0] -fastafile ciona_masked_transcripts.fasta -reldist 0 | wc -l";
    my $msg2 = `$todo`;  $msg2 =~ s/ //g; $msg2 =~ s/\n//g;

    my $todo = "/home/olly/PROGRAMS/GENREGEXP/genregexp -re $r->[0] -fastafile ciona.fasta -reldist 0 | wc -l";
    my $msg1 = `$todo`;  $msg1 =~ s/ //g; $msg1 =~ s/\n//g;

    my $ratio = $msg2 / $msg1;

    print "$r->[0]\t$msg2\t$msg1\t$ratio\n";
}
