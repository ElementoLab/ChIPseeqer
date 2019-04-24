use lib "$ENV{HOME}/PERL_MODULES";
use Repeats;
use Sets;
use Fasta;
use strict;

use Getopt::Long;

my $fastafile = undef;
my $extend    = 0;
my $output    = 0;

GetOptions("fastafile=s" => \$fastafile,
	   "output=s"    => \$output,
           "extend=s"    => \$extend);


my $fa = Fasta->new;
$fa->setFile($fastafile);

while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  
  my $tr = Repeats->new;
  $tr->setSeq($s);
  $tr->setProgram('TRF');
  $tr->setTRFminCopyNumber(5);
  $tr->setTRFmaxUnitSize(3);
  $tr->process();
  my $a_ref_tandems = $tr->getResults;
  $tr->dispose;
  
  my $tgtx = $extend;
  my $tgty = length($s) - $extend;

  if (@$a_ref_tandems > 0) {
    
    my $ov = 0;
    my $a_new_tandems = [];

    foreach my $r (@$a_ref_tandems) {
      #print "$r->[2]\n"; 
      if (Sets::sequencesOverlap($tgtx, $tgty, $r->[0], $r->[1])) {
	push @$a_new_tandems, $r;
      }
    }
    
    if (@$a_new_tandems > 0) {
      $s = lc($s);
      foreach my $r (@$a_new_tandems) {
	substr($s, $r->[0], $r->[1] - $r->[0] + 1) = uc(substr($s, $r->[0], $r->[1] - $r->[0] + 1));
      }
      
      if ($output == 0) {

	print ">$n\n";
	print "$s\n";

      } elsif ($output == 1) {
	
	my @b = split /\-/, $n;
	#print ">$n\n";
	print join("\t", @b) . "\t$s\n";
	
	
      }
    }
    
    #my $seq_masked = Sets::maskExons($s, $a_ref_tandems, 'X');
    
    
    #print "$seq_masked\n";
    #print "\n";
    #<STDIN>;
  }
}

