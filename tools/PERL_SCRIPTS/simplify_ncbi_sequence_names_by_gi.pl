#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    

    
    my ($gi) = $name =~ /^gi\|(\d+?)\|/;
    #$name =~ s/^.+\|//;



    print ">$gi\n$seq\n\n";
    
}
