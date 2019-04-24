#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use IO::File;

my $dir = $ARGV[1];
if (($dir ne "") && (! -e $dir)) {
  mkdir $dir;
}

open IN, $ARGV[0];
my $l = <IN>;
close IN;
my @a = split /\t/, $l, -1;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref   = $ta->getArray();
my $a_ref_1 = shift @$a_ref;
my $n1      = shift @$a_ref_1;
my $n       = scalar( @$a_ref_1 );
my @FH      = ();

foreach my $f (@$a_ref_1) {
  my $ff = $f;
  if (-e $f) {
    $ff = "$f-1";
  } 
  $ff =~ s/ /\_/g;
  $ff =~ s/[\)\(]//g;
  $ff =~ s/\|//g;
  $ff =~ s/\:/\_/g;
  $ff =~ s/\//\_/g;
  $ff =~ s/\*/\_/g;
  $ff =~ s/\"//g;
  $ff .= ".txt";
  if ($dir ne "") {
    $ff = "$dir/$ff";
    #$ff =~ s/\_.+$//;
  }

  

  my $fh = new IO::File $ff, "w";
  die "cannot create fh for $ff\n" if (!defined($fh)); 
  #print $fh "$n1\t$f\n";

  print $fh "GENE\tEXP\n";

  push @FH, $fh;
}


foreach my $r (@$a_ref) {
  my $m1 = shift @$r;
  for ($i=0; $i<$n; $i++) {
    my $fh = $FH[$i];
    print $fh "$m1\t$r->[$i]\n" if (($r->[$i] ne "") && (uc($r->[$i]) ne 'NAN'));
  }
}

foreach my $f (@FH) {
  undef $f;
}
  



