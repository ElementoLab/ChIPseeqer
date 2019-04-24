BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Fasta;
use Alignment;
use Sets;
use strict;

my $al = Alignment->new;
$al->setComp('RNA');

# read miRNA
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my $a_ref = $fa->nextSeq();
my ($n, $s) = @$a_ref;

$s = uc($s); 
$s =~ s/U/T/g;
my $mirna = $s; 

my @a = split //, $s;

$s = join("", reverse(@a));
 

my $ta = Table->new;

# read cluster membership

# go thru parition
$ta->loadFile($ARGV[2]);
my $h_ref = $ta->getIndexKV(0,1);

# go thru table
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {
  
  my $seq = uc($r->[3]);
  
  next if ($seq =~ /N/);
  next if (!defined($h_ref->{ $r->[0] }));

  my $a_ref_aln = $al->nw($seq, $s);

  print "$seq $mirna\n";


  my $sim = &rna_sim($a_ref_aln->[0], $a_ref_aln->[1]);
  push @{ $H{ $h_ref->{ $r->[0] } } }, $sim;

  if ($h_ref->{$r->[0]} == 4) {
    print "$r->[0]\t" . $h_ref->{$r->[0]} . "\n";
    #print "$seq\n   $s\n\n";
    print "$a_ref_aln->[0]\n$a_ref_aln->[1]\n\n";
    print "$sim";
    print "\n";
    print "-----------------------------------\n\n";
  }
}

foreach my $k (sort { $a <=> $b } keys(%H)) {
  #print "Cluster $k\n";
  my $a_ref = $H{$k};

  my $avg = Sets::average($a_ref);
  my $std = Sets::stddev($a_ref);

  print sprintf("Cluster %2d: avg match = %3.1f\% , std dev = %3.1f\%\n", $k, $avg*100, $std*100);

}


sub rna_sim {
  my ($s1, $s2) = @_;
  
  my @a1  = split //, $s1;
  my @a2  = split //, $s2;
  
  my $cnt = 0;
  my $cnt_t = 0;

  # advance up to last - on miRNA part
  my $i = 0;
  while ($a2[$i] eq '-') {
    #print "Advance $a2[$i]\n";
    $i ++;
  }

  for ($i=$i; $i<@a1; $i++) {

    my $n1 = $a1[$i];
    my $n2 = $a2[$i];
    if ((($n1 eq 'A') && ($n2 eq 'T')) || (($n1 eq 'T') && ($n2 eq 'A'))) {
      $cnt++;
    } elsif ((($n1 eq 'G') && ($n2 eq 'C')) || (($n1 eq 'C') && ($n2 eq 'G'))) {
      $cnt++;
    } elsif ((($n1 eq 'G') && ($n2 eq 'T')) || (($n1 eq 'T') && ($n2 eq 'G'))) {
      $cnt++;;
    } else {
      $cnt += 0;
    } 
    $cnt_t ++;
  }
  
  return $cnt / $cnt_t;
}
