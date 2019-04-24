my $pagedir ;
BEGIN{
    if ((!$ENV{PAGEDIR}) || ($ENV{PAGEDIR} eq '')) {
	$pagedir="./" ;
	print "The PAGEDIR environment variable is not set. It is set to default.\n";
    }
    else{
	$pagedir = $ENV{PAGEDIR};
    }
}

use lib "$pagedir/SCRIPTS";
use lib "$pagedir/SCRIPTS/PostScript-Simple-0.07/lib";

my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;

my $pvmatrixfile = undef ;
GetOptions ('pvmatrixfile=s'  => \$pvmatrixfile) ;

open LOG, "< $pvmatrixfile.log" or exit(0) ;

my %killed ;
while(<LOG>)
{
    s/\s+$// ;
    /^\[(.+)\].*\[(.+)\]$/ ;
    my ($c2, $c1) = ($1, $2) ;
    my ($p1, $p2) = split(/,/, $c1) ;
    $c1 = $p2.", ". $p1 ;
    my ($p1, $p2) = split(/,/, $c2) ;
    $c2 = $p2.", ". $p1 ;
    push(@{$killed{$c1}}, $c2) ;
}

open KILL, "> $pvmatrixfile.killed" or die "Couldn't open the .killed file" ;
my $dir = Sets::dirname($pvmatrixfile) ;
mkdir($dir."/SUMMARY") ;
foreach my $c1 (sort keys %killed){
    print KILL $c1, "\n" ;
    my $max = 0 ;
    my $cnt = @{$killed{$c1}} ;
    foreach my $c2 (@{$killed{$c1}}){
	print KILL "\t$c2\n" ;
	if (length($c2)>$max){
	    $max = length($c2) ;
	}
    }
    my $xbase        = 20 ;
    my $ybase        = 20 ;
    my $h            = 20 ;
    my $Font         = 17 ;
    my $font         = 14 ;
    my $space        = 50 ;
    my $ysize        = 2*$ybase + $h * $cnt ;
    my $xsize        = 2*$xbase + $max * $Font/2 + length($c1)* $font/2 + $space ;

    my $p = new PostScript::Simple(xsize     => $xsize,
                                   ysize     => $ysize,
                                   colour    => 1,
                                   eps       => 1,
                                   units     => "pt") ;
    
    $p->setlinewidth(1.0) ;
    $p->setcolour("black");
    $p->setfont("Garamond-bold", $Font);
    $p->text( { align => "right" }, $xbase + length($c1)* $Font/2, $ysize - ($ybase+$h * $cnt/2), $c1);
    my $y = 0 ;
    $p->setfont("Garamond", $font*1.0);

    foreach my $c2 (@{$killed{$c1}}){
	$p->setcolour("red");
	$p->line($xbase + length($c1)* $Font/2 + $font/2, $ysize - ($ybase + $h * $cnt/2 - $font/2), $xbase + length($c1)* $Font/2 + $space , $ysize - ($ybase + $y * $h - $font/2)) ;
	$p->setcolour("black");
	$p->text( { align => "left" }, $xbase + length($c1)* $Font/2 + $space + $font/2 , $ysize - ($ybase + $y * $h), $c2);
	$y++ ;
    }

    my ($dummy, $go) = split (/, /, $c1) ;
    my $outpdf = $dir."/SUMMARY/$go.pdf" ;
    $outpdf =~ s/\s/_/g ;
    my $outeps = $outpdf ;
    $outeps =~ s/pdf$/eps/ ;
    $p->output("$outeps") ;
    system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf") ;
}
