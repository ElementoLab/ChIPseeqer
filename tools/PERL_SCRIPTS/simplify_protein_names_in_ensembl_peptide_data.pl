#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    

    
    $name =~ s/[\t \|\(].+$//;
    #$name =~ s/^.+\|//;



    print ">$name\n$seq\n\n";
    
}
