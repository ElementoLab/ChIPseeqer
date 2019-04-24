#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
my $r     = shift @$a_ref;
#print     Sets::jointab($r);

print "$r->[0]\tEXP\n";

my $bkgkey = undef;
if ($ARGV[1] ne "") {
  $bkgkey = $ARGV[1];
}

my %H   = ();


foreach my $r (@$a_ref) {
  my @a = @$r;
  my $n = shift @a;
  my $s = join("/", @a);
  if (defined($bkgkey) && ($s eq $bkgkey)) {
  } else {
    $H{$s} = 1;
  }
}

my @sk = sort(keys(%H));
my $idx = 0;
my %NA = ();
if (defined($bkgkey)) {
  $NA{$bkgkey} = 0;
  $idx++;
}
foreach my $s (@sk) {
  $NA{$s} = $idx++;
}

foreach my $r (@$a_ref) {
  my $n = shift @$r;
  my $s = join("/", @$r);
  print "$n\t$NA{$s}\n";
}




foreach my $r (sort(keys(%NA))) {
  print STDERR "$r\t=> $NA{$r}\n";
}
