#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Superfamily;




print $ARGV[0] . "\n";

my $s = Superfamily->new;

my $a_ref = $s->getDomains($ARGV[0]);

my $a_new = $s->sortDomainsByRegion($a_ref);


foreach my $d (@$a_new) {
    print ">>" . $ARGV[0]; print "\t";
    print $d->{REGION};
    print "\t";
    print $d->{SFAMID};
    print "\t";
    print $d->{TEXT};
    print "\n";
}

