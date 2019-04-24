#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

shift @$a_ref;
my $cnt = 1;
my %C  = ();
foreach my $r (@$a_ref) {
  my $t = shift @$r;
  my $x = join(" ", @$r);
  $C{$x} ++;
}



my @o = (-1, 0, 1);

sub cmp_t {
  my ($t, $t1, $t2) = @_;
  if ($t < 0) {
    return "$t1 < $t2";
  } elsif ($t == 0) {
    return "$t1 ~ $t2";
  } elsif ($t > 0) {
    return "$t1 > $t2";
  }
}

my $cnt = 1;

my %NA = ();
my %N = ();
foreach my $r1 (@o) {
  my $t1 = &cmp_t($r1, "NS", "S");
  foreach my $r2 (@o) {
    my $t2 = &cmp_t($r2, "NS", "COPD");
    foreach my $r3 (@o) {
      my $t3 = &cmp_t($r3, "S", "COPD");
      my $k = "$r1 $r2 $r3";



      if ($C{$k} > 20) {
	print STDERR "$cnt => $k\n";
	$N{$k} = $cnt++;
	$NA{$k} = "$t1    $t2    $t3";
      }
    }
  }
}

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

shift @$a_ref;
my $cnt = 1;
my %C  = ();
foreach my $r (@$a_ref) {
  my $t = shift @$r;
  my $x = join(" ", @$r);

  if (defined($N{$x})) {
    print "$t\t$N{$x}\t$NA{$x}\n";    
  }
}


#foreach my $k (sort(keys(%H))) {
#  print STDERR "$H{$k} => $k \t($C{$k} genes)\n";
#}

