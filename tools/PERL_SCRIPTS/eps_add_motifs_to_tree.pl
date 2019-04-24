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


my $a_ref = Sets::getFiles("$ARGV[1]/*.eps");

foreach my $r (@$a_ref) {
  $r = Sets::filename($r);
  #my $c = $r;
  
}

foreach my $m (@mot) {
  
  my $them = $m->[0];
  my @a = split /\ /, $them;
  my $g = undef;
  my $d = undef;
  if (@a == 3) {
    $d = pop @a;
  }
  if ((@a == 2) && (($a[1] eq "A") || ($a[1] eq "B") || ($a[1] eq "C"))) {
    $d = pop @a;
  }
  $g = join("_", @a);

  my @matches = ();
  
  # go thru motifs
  foreach my $r (@$a_ref) {
    my $mog = $r;
    $mog =~ s/\_pwm\.txt.*$//;
    $mog =~ s/^.+?\_//;
    #$m =~ s/ext$//;
    #$m =~ s/D[123]$//;
    #$m =~ s/DsL//;
    $mog =~ s/MAL8P1153/MAL8P1\.153/;

    #my $mog = Sets::Pf_motifname2genename($r);
    print "Compare $mog with $g\n";
    if ($mog =~ /$g/) {
      
      print "ok, $mog contains $g\n";

      if (($mog =~ /DsL/) || !defined($d))  {
	push @matches, $r;
      } elsif (($mog =~ /D1/) && ($d eq "A")) {
	push @matches, $r;
      } elsif (($mog =~ /D2/) && ($d eq "B")) {
	push @matches, $r;
      } elsif (($mog =~ /D3/) && ($d eq "C")) {
	push @matches, $r;
      }
    }
  }
  
  my $pos = 0;
  foreach my $r (@matches) {
    my $e = new PostScript::Simple::EPS(file => "$ARGV[1]/$r");
    $e->scale(0.25);
    
    # add eps to the current page
    my $add = 0;
    if ($r =~ /GGGTGCACC/) {
      $add = 20;
    }	
    $p->importeps($e, $m->[1] + length($m->[0])*4 + $pos*40 + $add, $m->[2]-10);
    $pos++;
  }	


}

$p->output("file.ps");
system("ps2pdf file.ps");


exit(0);

my $cnt = 0;
foreach my $r (@$a_ref) {
  die if (! -e $r);
  my $e = new PostScript::Simple::EPS(file => $r);
  $e->scale(0.25);
  
  # add eps to the current page
  $p->importeps($e, $mot[$cnt]->[1] + length($mot[$cnt]->[0])*4, $mot[$cnt]->[2]-10);
  $cnt ++;
}
