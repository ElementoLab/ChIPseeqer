#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --f1=FILE --f2=FILE\n";
}

my $files = undef;
my $type  = "Alu";
my $fam   = 0;

GetOptions("files=s" => \$files,
	   "fam=s"   => \$fam,
	   "type=s"  => \$type);

my $a_ref = Sets::getFiles($files);


foreach my $f (@$a_ref) {
  my $todo = "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/CalcRepeatFamilyAbundanceInLibrary.pl $f";
  if ($fam == 1) {
    $todo .= " 1 ";
  }
  my $txt = `$todo`;
  my @a   = split /\n/, $txt;
  my %H   = ();
  foreach my $r (@a) {
    my @b = split /\t/, $r;
    $H{$b[0]} = $b[1];
  }
  my $ff = $f;
  $ff =~ s/\_combined\_mastertable.+$//;
  $ff =~ s/\_master\_table.+$//;

  print "$ff\t" . $H{$type} . "\n";
}
