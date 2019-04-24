#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use PostScript::Simple;
use Sets;
use strict;

# find text positions;

open IN, $ARGV[0];
my @lines = <IN>;

my @mot = ();
for (my $i=0; $i<@lines; $i++) {
  my $l = $lines[$i];
  if ($l =~ /^\((.+?)\) show/) {
    my $txt = $1;
    my $lin = $lines[$i-2];
    my @a   = split /\ /, $lin;
    print "$txt\t$a[0]\t$a[1]\n";

    # text, x, y
    my @a_tmp = ($txt, $a[0], $a[1]);
    push @mot, \@a_tmp;
  }
}
close IN;


# create a new PostScript object
my $p = new PostScript::Simple(papersize => "A4",
			    colour    => 1,
			    units     => "pt");

# create a new page
#$p->newpage;

# create an eps object
my $e = new PostScript::Simple::EPS(file => "$ARGV[0]");
#$e->rotate(90);
#$e->scale(0.5);

# add eps to the current page
$p->importeps($e, 0, 0);



foreach my $m (@mot) {

  my $add = 0;
  if ($m->[0] =~ /T1/) {
    $p->setcolour("red");
  } elsif ($m->[0] =~ /T2/) {
    $p->setcolour("green");
  } else {
    $p->setcolour("black");
  }
  $p->circle({filled => 1}, $m->[1] + length($m->[0])*4 + $add, $m->[2], 5);

}

$p->output("file.ps");
system("ps2pdf file.ps");


exit(0);
