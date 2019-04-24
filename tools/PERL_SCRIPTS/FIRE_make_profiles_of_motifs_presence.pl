#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Getopt::Long;
use Table;
use strict;

my $expfile = undef;
my $bin     = 1;

GetOptions("expfile=s" => \$expfile,
           "bin=s"     => \$bin);


my $file      = Sets::filename($expfile);
my $dir       = "$expfile\_FIRE/DNA/";
my $profile   = "$dir/$file.profiles";
my $sumfile   = "$dir/$file.summary";

my $ta = Table->new;

# load summary
$ta->loadFile($sumfile);
my $h_ref_en = $ta->getIndexKV(0,12);

# load profiles
$ta->loadFile($profile);
my %MOTIFS = ();
my %MATCHE = ();
my $a_ref_p = $ta->getArray();

foreach my $l (@$a_ref_p) {
  next if ($h_ref_en->{$l->[0]} != $bin);
  $MOTIFS{$l->[0]} = 1;
  $MATCHE{$l->[1]}{$l->[0]} ++;
}

$ta->loadFile($expfile);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;


my @mot = keys(%MOTIFS);
print "ID\t" . join("\t", @mot) . "\n";
foreach my $r (@$a_ref) {
  if ($r->[1] == $bin) {
    print "$r->[0]";
    foreach my $m (@mot) {
      my $c = $MATCHE{$r->[0]}{$m} * 1.0;
      print "\t$c";
    }
    print "\n";
  }
}

