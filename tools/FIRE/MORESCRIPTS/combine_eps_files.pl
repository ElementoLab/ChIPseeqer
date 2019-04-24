BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib qw(PostScript-Simple-0.07/lib);



use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;


my $eps1 = $ARGV[0];
my $eps2 = $ARGV[1];

my $e1  = new PostScript::Simple::EPS(file => $eps1);
my $eh1 = $e1->height;
my $ew1 = $e1->width;

my $e2  = new PostScript::Simple::EPS(file => $eps2);
my $eh2 = $e2->height;
my $ew2 = $e2->width;

my $ysize = $eh1+$eh2-700;

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


$p->_add_eps($e2, 0, $ysize-$eh1-$eh2+360);
 


$p->output("out.eps");
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None out.eps out.pdf");




