#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();
my $h_ref = {};
foreach my $r (@$a_ref) {
  $r->[0] =~ s/\ //g;
  $h_ref->{$r->[0]} = $r->[1];
}


my $w = undef;
my $neww = 26;

my $patient = undef;
if ($ARGV[2] ne "") {
  $patient = $ARGV[2];
}

open IN, $ARGV[0];
my $l = <IN>;
print $l;
while (my $l = <IN>) {
  chomp $l;
  if ($l eq "") {
    print "$l\n";
  } else {	

    if (!defined($w)) {
      my ($nt) = $l =~ /^(.+?\ +)/;
      $w = length($nt);
    }

    my ($n1, $s) = $l =~ /^(.{$w})(.+)$/;
    
    $n1 =~ s/\ //g;
    
    if (defined($h_ref->{$n1})) {
      $n1 = "$n1-" . $h_ref->{$n1};

      if (defined($patient)) {
	$n1 = "$patient-$n1";
      }

    } else {
      print STDERR "Cannot find name for '$n1'\n";
      $n1 = "$n1";
    }	

    print $n1 . (" " x ($neww-length($n1))) . "$s\n";
  }
  
}
close IN;
