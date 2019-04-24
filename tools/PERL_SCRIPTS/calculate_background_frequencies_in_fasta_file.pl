#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);


my %H = ();
my $cnt = 0;
my %N = ();
while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    #next if (length($seq) < 50);
    
    my $seq = uc($seq);

    my @a = split //, $seq;

    foreach my $n (@a) {
        if ($n =~ /[ACTG]/i) {
	  $H{$n}++; $cnt++;
        } else {
          $N{$n}++;
        }
    }
    
    

    #print ">$name\n$seq\n\n";
    
}


foreach my $k (keys(%H)) {
    print sprintf("$k\t%4.3f\t%d\t%d\n", 1.0*$H{$k}/$cnt, $H{$k}, $cnt);
}
print "G+C:\n";
my $gc = ($H{"C"}+$H{"G"})/$cnt;
print sprintf("%4.3f\n", $gc);
print "Other characters:\n";
foreach my $k (keys(%N)) {
   print "$k\t$N{$k}\n";
}
