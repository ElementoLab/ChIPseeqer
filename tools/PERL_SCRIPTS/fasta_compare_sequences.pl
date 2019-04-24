#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;

my %SEQ = ();
my $fa1 = Fasta->new;
$fa1->setFile($ARGV[0]);

while (my $a_ref = $fa1->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $n =~ s/\r//g;
    $SEQ{$n} = $s;
    
}



my $fa2 = Fasta->new;
$fa2->setFile($ARGV[1]);

while (my $a_ref = $fa2->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $n =~ s/\r//g;

    if (!defined($SEQ{$n})) {
      print "Cannot find $n\n";
    } else {
      
      if ($s eq $SEQ{$n}) {
	print "Sequences $n are the same in two files\n";
      } else {
	print "Sequences $n are different.\n";
	print "$SEQ{$n}\n";
	print "$s\n";
	
      }
      
    }


}
