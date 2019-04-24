#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}/SCRIPTS/PostScript-Simple-0.07/lib";
use lib "$ENV{HOME}/PERL_MODULES";

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;


my $eps1 = $ARGV[0];
my $eps2 = $ARGV[1];

my $up2  = 0;
my $outfile = "out.eps";
my $label = undef;
GetOptions("up2=s"     => \$up2,
	   "label=s"   => \$label,
           "outfile=s" => \$outfile);



my $e1  = new PostScript::Simple::EPS(file => $eps1);
my $eh1 = $e1->height;
my $ew1 = $e1->width;

my $e2  = new PostScript::Simple::EPS(file => $eps2);
my $eh2 = $e2->height;
my $ew2 = $e2->width;

my $ysize = $eh1+$eh2-$up2;


my $p = new PostScript::Simple(xsize => $ew1, #papersize => "A4",
			       ysize => $ysize, #papersize => "A4",
			       #papersize => "A4",
			       colour    => 1,
			       #landscape    => 1,
			       #coordorigin => "LeftTop",
			       #direction   => "RightDown",
			       eps       => 1,
			       units     => "pt");



$p->_add_eps($e1, 0, $ysize-$eh1);


$p->_add_eps($e2, 0, $ysize-($eh1+$eh2-$up2));

if (defined($label)) {
  $p->setcolour("black");
  $p->setfont("Arial", 15);
  $p->text({ align => 'left', rotate => 0 }, 10, $ysize-20, $label);
} 


$p->output($outfile);
my $outpdf = $outfile;
$outpdf =~ s/\.eps$/\.pdf/;
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outfile $outpdf");

if (-e $outpdf) {
  print "Created $outpdf\n";
  system("evince $outpdf");
}



