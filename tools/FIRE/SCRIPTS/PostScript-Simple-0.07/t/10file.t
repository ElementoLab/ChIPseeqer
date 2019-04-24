#!/usr/bin/perl -w

use strict;
use lib qw(./lib ../lib t/lib);
use Test::Simple tests => 7;
#use Data::Dumper;
use PostScript::Simple;

my $f = 'xtest-a.ps';
my $s = new PostScript::Simple(xsize => 50, ysize => 200);

$s->box(10, 10, 40, 190);
$s->output( $f );

#print STDERR Dumper $s;

# check object
ok( $s->{usedbox} == 1 );
ok( $s->{psfunctions} =~ m|/u | );
ok( index( $s->{pspages}, q[10 ux 10 uy 40 ux 190 uy box stroke]) > -1 );

# check output
ok( -e $f );
open( CHK, $f ) or die("Can't open the file $f: $!");
$/ = undef;
my $file = <CHK>;
close CHK;

ok( index( $file, '%!PS-Adobe-3.0 EPSF-1.2' ) == 0 );
ok( index( $file, '%%EOF' ) == (length( $file ) - 6) );
ok( index( $file, '10 ux 10 uy 40 ux 190 uy box stroke' ) > 0 );
