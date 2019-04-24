#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

my $min = (defined($ARGV[1])?$ARGV[1]:50);

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    next if (length($seq) < $min);
    
    $name =~ s/\ .+$//;
    
    print ">$name\n$seq\n\n";
    
}
