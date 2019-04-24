#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);


while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    my @a = split /[\|\ ]/, $name;

    
    
    print ">$a[1]\n$seq\n\n";
    
}
