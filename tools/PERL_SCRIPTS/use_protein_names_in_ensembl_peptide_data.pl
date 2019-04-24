#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

my $pat = undef;
if ($ARGV[1]) {
  $pat = $ARGV[1];
} else {
  $pat = "ENSP";
}

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    

    
    $name =~ /\|($pat.+?)\ /;
    
    print ">$1\n$seq\n\n";
    
}
