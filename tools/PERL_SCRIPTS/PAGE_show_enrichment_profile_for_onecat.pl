#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Getopt::Long;



# load onecat


# traverse expfile, via windows 
if (@ARGV == 0) {
  die "Args : --expfile --onecatfile\n";
}

GetOptions("expfile=s"    => \$expfile,
	   "onecatfile=s" => \$onecatfile);

# load expfile
my $ta = Table->new;
$ta->loadFile($expfile);

my $a_ref = $ta->getArray();
shift @$a_ref;

# sort by col 1
$ta->sortbycol(1);

# load onecatfile
$ta->loadFile($onecatfile);
my $h_ref = $ta->getIndexKV(0,1);

# need to remove genes for which we don't have annot
my @new_genes = ();
foreach my $r (@$a_ref) {
  push @new_genes, $r if (defined($h_ref->{$r->[0]}));
}

$a_ref = \@new_genes;

my $n = @$a_ref;
my $w = 100;

for (my $i=0; $i<$n-$w; $i++) {
  #print STDERR "$a_ref->[$i]->[1]\n";
  my $cnt = 0;
  for (my $j=$i; $j<$i+$w; $j++) {
    if ($h_ref->{$a_ref->[$j]->[0]} == 1) {
      $cnt ++;
    }
  }

  my $f = $cnt/$w;
  print "$f\n";

}

# count
my $frac = 0;
for (my $i=0; $i<$n; $i++) {

  if ($h_ref->{$a_ref->[$i]->[0]} == 1) {
    $frac ++;
  }
}

$frac /= $n;

print STDERR "FRAC=$frac\n";






