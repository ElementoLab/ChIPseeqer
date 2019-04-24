#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref1 = $ta->getColumn($ARGV[0]);
my $a_ref3 = $ta->getArray(); 

$ta->loadFile($ARGV[2]);
my $a_ref2 = $ta->getColumn($ARGV[0]);

my $a_union = Sets::getOverlapSet($a_ref1, $a_ref2);

my %H = ();
foreach my $h (@$a_union) {
    $H{ $h } = 1;
}

foreach my $r (@$a_ref3) {

    if (defined($H{ $r->[ $ARGV[0] ] })) {

	print join("\t", @$r); print "\n";

    } 

}
