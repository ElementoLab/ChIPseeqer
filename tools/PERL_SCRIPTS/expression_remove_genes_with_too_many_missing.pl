#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
if (@ARGV == 0) {
  die "Args: table maxnummissing\n";
}
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";
foreach my $r (@$a_ref) {
  my $m = shift @$r;
  my $cnt = 0;
  my $n = scalar(@$r);
  foreach my $s (@$r) {
    $s =~ s/\ //g;
    if (($s eq "") || ($s eq "NA")) {
      $cnt ++;
    }
  }

  if ($cnt >= $ARGV[1]) {
    print STDERR "Skipping $m\n";
  } else {
    print "$m\t" . join("\t", @$r) . "\n";
  }
  

}

