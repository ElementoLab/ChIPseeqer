BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Fasta;
use Sets;
use Vienna;
use strict;


my $vi = Vienna->new;

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


  next if (length($seq) != 25);

  if ($ARGV[3]) {
    # strip off nt around seed match
    my ($p1, $p2) = $seq =~ /^(.+)(.{8})$/;
    my @a_p2 = split //, $p2;  shift @a_p2; pop @a_p2; $p2 = join('', @a_p2);
    $seq = "$p1$p2";
  }
  

  my $thismirna = $mirna;

  if ($ARGV[3]) {
    # strip off nt around mirna seed
    #print "BEF: $thismirna\n";
    my ($p1, $p2) = $thismirna =~ /^(.{8})(.+)$/;
    my @a_p1 = split //, $p1;  shift @a_p1; pop @a_p1; $p1 = join('', @a_p1);
    $thismirna = "$p1$p2";
    #print "AFT: $thismirna\n";
  }
  





  $vi->setCofoldSequences($seq, $thismirna);

  my $a_ref_folds = $vi->fold();


  my $first_fold  = shift @$a_ref_folds;
  
  #next if (length($first_fold->{SEQUENCE}) != 48);

  #

  if ($ARGV[3]) {
    next if ($first_fold->{FOLD} !~ /\(\(\(\(\(\(\&\)\)\)\)\)\)/);
  } else {
    next if ($first_fold->{FOLD} !~ /\(\(\(\(\(\([\.\(]\&[\.\)]\)\)\)\)\)\)/);
  }

  #print "E=$first_fold->{ENERGY}\n";
  
  my $sim         = $first_fold->{ENERGY};

  push @{ $H{ $h_ref->{ $r->[0] } } }, $sim;

  #if ($h_ref->{$r->[0]} == 4) {
  #  print $vi->getRawOutput();
  #  print "-------------------------------\n";
  #}

  
}

# output to disk


foreach my $k (sort { $a <=> $b } keys(%H)) {
  #print "Cluster $k\n";
  my $a_ref = $H{$k};

  Sets::writeSet($a_ref, "folds_C$k.txt");


  my $avg = Sets::average($a_ref);
  my $std = Sets::stddev($a_ref);

  print sprintf("Cluster %2d: avg match = %3.1f\% , std dev = %3.1f\%\n", $k, $avg, $std);

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
