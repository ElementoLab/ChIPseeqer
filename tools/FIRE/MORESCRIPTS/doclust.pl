use AggloClust;
use Table;
use Getopt::Long;
use strict;

my $cor     = undef;
my $dat     = undef;
my $fulldat = undef;
my $mind    = 0.75;
GetOptions ('cor=s'     => \$cor,
	    'mind=s'    => \$mind,
	    'dat=s'     => \$dat,
	    'fulldat=s' => \$fulldat);


my $ta = Table->new;
$ta->loadFile($cor);
my $a_ref = $ta->getArray();
shift @$a_ref;
foreach my $r (@$a_ref) {
  shift @$r;
  foreach my $c (@$r) {
    $c = 1 - $c;
  }
}





my $ag = AggloClust->new;
$ag->setDistanceMatrix($a_ref);
$ag->setMinD($mind);
#$ag->setMaxNbClusters(10);

my $a_ref_c = $ag->agglomerate_using_avg_linkage();


$ta->loadFile($dat);
my $a_ref_d = $ta->getArray();


foreach my $r (@$a_ref_c) {
  print join(" ", @$r); print "\n";
}


my @GENES = ();
open IN, $fulldat;
my $l = <IN>;
while (my $l = <IN>) {
  my @a = split /\t/, $l;
  push @GENES, $a[0];
}

my $ff = $dat; $ff =~ s/\.txt/\.CDT/;
die "toto" if ($ff eq $dat);
open OUT, ">$ff";
my $r = shift @$a_ref_d;
my $n = shift @$r;
my %H = ();
print OUT "I\tNAME\tGWEIGHT\t" . join("\t", @$r) . "\n";
print OUT "EWEIGHT\t\t" . ("\t1.000" x @$r) . "\n";
my $cnt = 0;
foreach my $r (@$a_ref_c) {
  foreach my $s (@$r) {
    my $n = shift @{ $a_ref_d->[$s] };
    $H{ $n } = $cnt;
    print OUT "$n\t$n\t1.000\t" . join("\t", @{ $a_ref_d->[$s] }) . "\n";
  }
  $cnt ++;
}
close OUT;


# create partition
open OUT, ">part.txt";
print OUT "i\ti\n";
foreach my $g (@GENES) {
  print OUT "$g\t";
  if (defined($H{$g})) {
    print OUT "$H{$g}\n";
  } else {
    print OUT "$cnt\n";
  }
}
close OUT;
