#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    

    
    $name =~ s/\|.+$//;
    $name =~ s/^.+\|//;
    $name =~ s/\.1$//;


    print ">$name\n$seq\n\n";
    
}
