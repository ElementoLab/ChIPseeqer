#!/usr/bin/perl -w

use strict;
use lib qw(./lib ../lib t/lib);
use Test::Simple tests => 11;
#use Data::Dumper;
use PostScript::Simple;

my $f = "xtest-b.ps";
my $t = new PostScript::Simple(landscape => 0,
            eps => 0,
            papersize => "a4",
            colour => 1,
            clip => 0,
            units => "mm");

ok( $t );

$t->newpage(-1);

$t->line(10,10, 10,50);
$t->setlinewidth(8);
$t->line(90,10, 90,50);
$t->linextend(40,90);
$t->setcolour("brightred");
$t->circle({filled=>1}, 40, 90, 30);
$t->setcolour("darkgreen");
$t->setlinewidth(0.1);
for (my $i = 0; $i < 360; $i += 20) {
  $t->polygon({offset=>[0,0], rotate=>[$i,70,90], filled=>0}, 40,90, 69,92, 75,84);
}

$t->setlinewidth("thin");
$t->setcolour("darkgreen");
$t->box(20, 10, 80, 20);
$t->setcolour("grey30");
$t->box({filled=>1}, 20, 30, 80, 40);
$t->setcolour("grey10");
$t->setfont("Bookman", 12);
$t->text(5,5, "Matthew");

$t->newpage;
$t->line((10, 20), (30, 40));
$t->linextend(60, 50);

$t->line(10,12, 20,12);
$t->polygon(10,10, 20,10);

$t->setcolour("grey90");
$t->polygon({offset=>[5,5], filled=>1}, 10,10, 15,20, 25,20, 30,10, 15,15, 10,10, 0);
$t->setcolour("black");
$t->polygon({offset=>[10,10], rotate=>[45,20,20]}, 10,10, 15,20, 25,20, 30,10, 15,15, 10,10, 1);

$t->line((0, 100), (100, 0), (255, 0, 0));

$t->newpage(30);

for (my $i = 12; $i < 80; $i += 2) {
  $t->setcolour($i*3, 0, 0);
  $t->box({filled=>1}, $i - 2, 10, $i, 40);
}

$t->line((40, 30), (30, 10));
$t->linextend(60, 0);
$t->line((0, 100), (100, 0),(0, 255, 0));

$t->output( $f );
#$t->output( "x" );

ok( -e $f );

open( FILE, $f ) or die("Can't open $f: $!");
$/ = undef;
my $lines = <FILE>;
close FILE;

ok( $lines =~ m/%%LanguageLevel: 1/s );
ok( $lines =~ m/%%DocumentMedia: A4 595.27559 841.88976 0 \( \) \( \)/s );
ok( $lines =~ m/%%Orientation: Portrait/s );
ok( $lines =~ m/%%Pages: 3/s );

ok( index($lines, "%!PS-Adobe-3.0\n") == 0 );
my ( $prolog ) = ( $lines =~ m/%%BeginResource: PostScript::Simple\n(.*)%%EndResource/s );
ok( $prolog );
ok( $prolog eq PROLOG());

my ( $body ) = ( $lines =~ m/%%EndProlog\n(.*)%%EOF/s );
ok( $body );
ok( $body eq BODY());

#print ">>>$body<<<<<<\n";

### Subs

sub PROLOG {
	return q[/ux {72 mul 25.4 div} def
/uy {72 mul 25.4 div} def
/u {72 mul 25.4 div} def
/STARTDIFFENC { mark } bind def
/ENDDIFFENC { 

% /NewEnc BaseEnc STARTDIFFENC number or glyphname ... ENDDIFFENC -
	counttomark 2 add -1 roll 256 array copy
	/TempEncode exch def
	
	% pointer for sequential encodings
	/EncodePointer 0 def
	{
		% Get the bottom object
		counttomark -1 roll
		% Is it a mark?
		dup type dup /marktype eq {
			% End of encoding
			pop pop exit
		} {
			/nametype eq {
			% Insert the name at EncodePointer 

			% and increment the pointer.
			TempEncode EncodePointer 3 -1 roll put
			/EncodePointer EncodePointer 1 add def
			} {
			% Set the EncodePointer to the number
			/EncodePointer exch def
			} ifelse
		} ifelse
	} loop	

	TempEncode def
} bind def

% Define ISO Latin1 encoding if it doesnt exist
/ISOLatin1Encoding where {
%	(ISOLatin1 exists!) =
	pop
} {
	(ISOLatin1 does not exist, creating...) =
	/ISOLatin1Encoding StandardEncoding STARTDIFFENC
		144 /dotlessi /grave /acute /circumflex /tilde 
		/macron /breve /dotaccent /dieresis /.notdef /ring 
		/cedilla /.notdef /hungarumlaut /ogonek /caron /space 
		/exclamdown /cent /sterling /currency /yen /brokenbar 
		/section /dieresis /copyright /ordfeminine 
		/guillemotleft /logicalnot /hyphen /registered 
		/macron /degree /plusminus /twosuperior 
		/threesuperior /acute /mu /paragraph /periodcentered 
		/cedilla /onesuperior /ordmasculine /guillemotright 
		/onequarter /onehalf /threequarters /questiondown 
		/Agrave /Aacute /Acircumflex /Atilde /Adieresis 
		/Aring /AE /Ccedilla /Egrave /Eacute /Ecircumflex 
		/Edieresis /Igrave /Iacute /Icircumflex /Idieresis 
		/Eth /Ntilde /Ograve /Oacute /Ocircumflex /Otilde 
		/Odieresis /multiply /Oslash /Ugrave /Uacute 
		/Ucircumflex /Udieresis /Yacute /Thorn /germandbls 
		/agrave /aacute /acircumflex /atilde /adieresis 
		/aring /ae /ccedilla /egrave /eacute /ecircumflex 
		/edieresis /igrave /iacute /icircumflex /idieresis 
		/eth /ntilde /ograve /oacute /ocircumflex /otilde 
		/odieresis /divide /oslash /ugrave /uacute 
		/ucircumflex /udieresis /yacute /thorn /ydieresis
	ENDDIFFENC
} ifelse

% Name: Re-encode Font
% Description: Creates a new font using the named encoding. 

/REENCODEFONT { % /Newfont NewEncoding /Oldfont
	findfont dup length 4 add dict
	begin
		{ % forall
			1 index /FID ne 
			2 index /UniqueID ne and
			2 index /XUID ne and
			{ def } { pop pop } ifelse
		} forall
		/Encoding exch def
		% defs for DPS
		/BitmapWidths false def
		/ExactSize 0 def
		/InBetweenSize 0 def
		/TransformedChar 0 def
		currentdict
	end
	definefont pop
} bind def

% Reencode the std fonts: 
/Courier-iso ISOLatin1Encoding /Courier REENCODEFONT
/Courier-Bold-iso ISOLatin1Encoding /Courier-Bold REENCODEFONT
/Courier-BoldOblique-iso ISOLatin1Encoding /Courier-BoldOblique REENCODEFONT
/Courier-Oblique-iso ISOLatin1Encoding /Courier-Oblique REENCODEFONT
/Helvetica-iso ISOLatin1Encoding /Helvetica REENCODEFONT
/Helvetica-Bold-iso ISOLatin1Encoding /Helvetica-Bold REENCODEFONT
/Helvetica-BoldOblique-iso ISOLatin1Encoding /Helvetica-BoldOblique REENCODEFONT
/Helvetica-Oblique-iso ISOLatin1Encoding /Helvetica-Oblique REENCODEFONT
/Times-Roman-iso ISOLatin1Encoding /Times-Roman REENCODEFONT
/Times-Bold-iso ISOLatin1Encoding /Times-Bold REENCODEFONT
/Times-BoldItalic-iso ISOLatin1Encoding /Times-BoldItalic REENCODEFONT
/Times-Italic-iso ISOLatin1Encoding /Times-Italic REENCODEFONT
/Symbol-iso ISOLatin1Encoding /Symbol REENCODEFONT
/circle {newpath 0 360 arc closepath} bind def
/rotabout {3 copy pop translate rotate exch 0 exch
sub exch 0 exch sub translate} def
/box {
  newpath 3 copy pop exch 4 copy pop pop
  8 copy pop pop pop pop exch pop exch
  3 copy pop pop exch moveto lineto
  lineto lineto pop pop pop pop closepath
} bind def
];
}

sub BODY {
	return q[%%BeginSetup
ll 2 ge { << /PageSize [ 595.27559 841.88976 ] /ImagingBBox null >> setpagedevice } if
%%EndSetup
%%Page: -1 1
%%BeginPageSetup
/pagelevel save def
%%EndPageSetup
newpath
10 ux 10 uy moveto
10 ux 50 uy lineto stroke
8 u setlinewidth
newpath
90 ux 10 uy moveto
90 ux 50 uy lineto
40 ux 90 uy lineto stroke
1 0 0 setrgbcolor
40 ux 90 uy 30 u circle fill
0 0.5 0 setrgbcolor
0.1 u setlinewidth
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
gsave 70 ux 90 uy 20 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 40 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 60 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 80 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 100 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 120 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 140 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 160 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 180 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 200 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 220 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 240 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 260 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 280 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 300 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 320 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
gsave 70 ux 90 uy 340 rotabout
newpath
40 ux 90 uy moveto
69 ux 92 uy lineto 75 ux 84 uy lineto stroke
grestore
0.4 setlinewidth
0 0.5 0 setrgbcolor
20 ux 10 uy 80 ux 20 uy box stroke
0.3 0.3 0.3 setrgbcolor
20 ux 30 uy 80 ux 40 uy box fill
0.1 0.1 0.1 setrgbcolor
/Bookman findfont 12 scalefont setfont
newpath
5 ux 5 uy moveto
(Matthew)   show stroke 
%%PageTrailer
pagelevel restore
showpage
%%Page: -2 2
%%BeginPageSetup
/pagelevel save def
%%EndPageSetup
newpath
10 ux 20 uy moveto
30 ux 40 uy lineto
60 ux 50 uy lineto stroke
newpath
10 ux 12 uy moveto
20 ux 12 uy lineto stroke
newpath
10 ux 10 uy moveto
20 ux 10 uy lineto stroke
0.9 0.9 0.9 setrgbcolor
gsave 5 ux 5 uy translate
newpath
10 ux 10 uy moveto
15 ux 20 uy lineto 25 ux 20 uy lineto 30 ux 10 uy lineto 15 ux 15 uy lineto 10 ux 10 uy lineto fill
grestore
0 0 0 setrgbcolor
gsave 10 ux 10 uy translate
20 ux 20 uy 45 rotabout
newpath
10 ux 10 uy moveto
15 ux 20 uy lineto 25 ux 20 uy lineto 30 ux 10 uy lineto 15 ux 15 uy lineto 10 ux 10 uy lineto stroke
grestore
1 0 0 setrgbcolor
newpath
0 ux 100 uy moveto
100 ux 0 uy lineto stroke
%%PageTrailer
pagelevel restore
showpage
%%Page: 30 3
%%BeginPageSetup
/pagelevel save def
%%EndPageSetup
0.141176470588235 0 0 setrgbcolor
10 ux 10 uy 12 ux 40 uy box fill
0.164705882352941 0 0 setrgbcolor
12 ux 10 uy 14 ux 40 uy box fill
0.188235294117647 0 0 setrgbcolor
14 ux 10 uy 16 ux 40 uy box fill
0.211764705882353 0 0 setrgbcolor
16 ux 10 uy 18 ux 40 uy box fill
0.235294117647059 0 0 setrgbcolor
18 ux 10 uy 20 ux 40 uy box fill
0.258823529411765 0 0 setrgbcolor
20 ux 10 uy 22 ux 40 uy box fill
0.282352941176471 0 0 setrgbcolor
22 ux 10 uy 24 ux 40 uy box fill
0.305882352941176 0 0 setrgbcolor
24 ux 10 uy 26 ux 40 uy box fill
0.329411764705882 0 0 setrgbcolor
26 ux 10 uy 28 ux 40 uy box fill
0.352941176470588 0 0 setrgbcolor
28 ux 10 uy 30 ux 40 uy box fill
0.376470588235294 0 0 setrgbcolor
30 ux 10 uy 32 ux 40 uy box fill
0.4 0 0 setrgbcolor
32 ux 10 uy 34 ux 40 uy box fill
0.423529411764706 0 0 setrgbcolor
34 ux 10 uy 36 ux 40 uy box fill
0.447058823529412 0 0 setrgbcolor
36 ux 10 uy 38 ux 40 uy box fill
0.470588235294118 0 0 setrgbcolor
38 ux 10 uy 40 ux 40 uy box fill
0.494117647058824 0 0 setrgbcolor
40 ux 10 uy 42 ux 40 uy box fill
0.517647058823529 0 0 setrgbcolor
42 ux 10 uy 44 ux 40 uy box fill
0.541176470588235 0 0 setrgbcolor
44 ux 10 uy 46 ux 40 uy box fill
0.564705882352941 0 0 setrgbcolor
46 ux 10 uy 48 ux 40 uy box fill
0.588235294117647 0 0 setrgbcolor
48 ux 10 uy 50 ux 40 uy box fill
0.611764705882353 0 0 setrgbcolor
50 ux 10 uy 52 ux 40 uy box fill
0.635294117647059 0 0 setrgbcolor
52 ux 10 uy 54 ux 40 uy box fill
0.658823529411765 0 0 setrgbcolor
54 ux 10 uy 56 ux 40 uy box fill
0.682352941176471 0 0 setrgbcolor
56 ux 10 uy 58 ux 40 uy box fill
0.705882352941177 0 0 setrgbcolor
58 ux 10 uy 60 ux 40 uy box fill
0.729411764705882 0 0 setrgbcolor
60 ux 10 uy 62 ux 40 uy box fill
0.752941176470588 0 0 setrgbcolor
62 ux 10 uy 64 ux 40 uy box fill
0.776470588235294 0 0 setrgbcolor
64 ux 10 uy 66 ux 40 uy box fill
0.8 0 0 setrgbcolor
66 ux 10 uy 68 ux 40 uy box fill
0.823529411764706 0 0 setrgbcolor
68 ux 10 uy 70 ux 40 uy box fill
0.847058823529412 0 0 setrgbcolor
70 ux 10 uy 72 ux 40 uy box fill
0.870588235294118 0 0 setrgbcolor
72 ux 10 uy 74 ux 40 uy box fill
0.894117647058824 0 0 setrgbcolor
74 ux 10 uy 76 ux 40 uy box fill
0.917647058823529 0 0 setrgbcolor
76 ux 10 uy 78 ux 40 uy box fill
newpath
40 ux 30 uy moveto
30 ux 10 uy lineto
60 ux 0 uy lineto stroke
0 1 0 setrgbcolor
newpath
0 ux 100 uy moveto
100 ux 0 uy lineto stroke
%%PageTrailer
pagelevel restore
showpage
];
}


