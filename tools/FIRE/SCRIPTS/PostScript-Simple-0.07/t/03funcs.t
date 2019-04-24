#!/usr/bin/perl
use strict;
use lib qw(./lib ../lib t/lib);
use Test::Simple tests => 44;
#use Data::Dumper;
use PostScript::Simple;

# huge workout of all methods, OK and error conditions

my $s = new PostScript::Simple(xsize => 350, ysize => 350, eps => 1, colour => 1);

ok( $s );
ok( ! $s->newpage );

eval { $s->output; };
ok( $@ );

ok( $s->setcolour('black') );
ok( $s->setcolour('BLACK') );
ok( ! $s->setcolour('Geddy lee') );
ok( ! $s->setcolour(120, 240) );
ok( $s->setcolour(120, 240, 0) );


ok( $s->setlinewidth(1) );
ok( ! $s->setlinewidth );


ok( $s->line(10,10, 10,20) );
ok( ! $s->line(10,10, 10,20, 50, 50) );
ok( ! $s->line(10,10, 10) );
ok( $s->line(10,10, 10,20, 50, 50, 50) );


ok( $s->linextend(100,100) );
ok( ! $s->linextend(100) );


ok( $s->polygon(10,10, 10,20, 110,10, 110,20) );
#ok( $s->polygon(10,10, 10,20, 110,10, 110,20, 1) );
ok( $s->polygon({rotate=>45,filled=>1}, 10,10, 10,20, 110,10, 110,20) );
ok( $s->polygon({rotate=>[45,20,20]}, 10,10, 10,20, 110,10, 110,20) );
ok( $s->polygon({offset=>[10,10]}, 10,10, 10,20, 110,10, 110,20) );
ok( ! $s->polygon(10,10, 10) );


ok( $s->circle( 120, 120, 30 ) );
ok( $s->circle( {filled=>1}, 120, 120, 30 ) );
ok( ! $s->circle( 120 ) );
ok( ! $s->circle );


ok( $s->box(210,210, 220,230) );
ok( $s->box( {filled=>1}, 215,215, 225,235) );
ok( ! $s->box(210,210, 220) );


ok( $s->setfont('Helvetica', 12) );
ok( ! $s->setfont('Helvetica') );


ok( $s->text( 10, 10, 'Hello World' ) );
ok( $s->text( {align=>"left"}, 10, 10, 'Hello World' ) );
ok( $s->text( {rotate=>56}, 10, 10, 'Hello World' ) );
ok( ! $s->text( 10, 10, 'Hello World', 'foo', 'wobble' ) );
ok( ! $s->text( 10, 10 ) );


ok( ! $s->curve(10,310, 10,320, 110,310, 110) );
ok( $s->curve(10,310, 10,320, 110,310, 110,320) );


ok( $s->curvextend(110,330, 210,330, 210,320) );
ok( ! $s->curvextend(110,330, 210,330, 210) );


ok( length($s->{'pspages'}) eq length(CANNED()) );
ok( $s->{'pspages'} eq CANNED() );

ok( length($s->{'psfunctions'}) eq length(FUNCS()) );
ok( $s->{'psfunctions'} eq FUNCS() );

ok( $s->output('x03.eps') );
unlink 'x03.eps';

#print Dumper $s;

###

sub FUNCS {
return '/ux {} def
/uy {} def
/u {} def
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
/rotabout {3 copy pop translate rotate exch 0 exch
sub exch 0 exch sub translate} def
/circle {newpath 0 360 arc closepath} bind def
/box {
  newpath 3 copy pop exch 4 copy pop pop
  8 copy pop pop pop pop exch pop exch
  3 copy pop pop exch moveto lineto
  lineto lineto pop pop pop pop closepath
} bind def
';
}

sub CANNED {
return '(error: Do not use newpage for eps files!
) print flush
0 0 0 setrgbcolor
0 0 0 setrgbcolor
(error: bad colour name \'geddy lee\'
) print flush
(error: setcolour given invalid arguments: 120, 240, undef
) print flush
0.470588235294118 0.941176470588235 0 setrgbcolor
1 u setlinewidth
(error: setlinewidth not given a width
) print flush
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto stroke
(error: wrong number of args for line
) print flush
(error: wrong number of args for line
) print flush
0.196078431372549 0.196078431372549 0.196078431372549 setrgbcolor
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto
100 ux 100 uy lineto stroke
(error: wrong number of args for linextend
) print flush
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto 110 ux 10 uy lineto 110 ux 20 uy lineto stroke
gsave 10 ux 10 uy 45 rotabout
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto 110 ux 10 uy lineto 110 ux 20 uy lineto fill
grestore
gsave 20 ux 20 uy 45 rotabout
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto 110 ux 10 uy lineto 110 ux 20 uy lineto stroke
grestore
gsave 10 ux 10 uy translate
newpath
10 ux 10 uy moveto
10 ux 20 uy lineto 110 ux 10 uy lineto 110 ux 20 uy lineto stroke
grestore
(error: bad polygon - not enough points
) print flush
120 ux 120 uy 30 u circle stroke
120 ux 120 uy 30 u circle fill
(error: circle: wrong number of arguments
) print flush
(error: circle: wrong number of arguments
) print flush
210 ux 210 uy 220 ux 230 uy box stroke
215 ux 215 uy 225 ux 235 uy box fill
(error: box: wrong number of arguments
) print flush
/Helvetica findfont 12 scalefont setfont
(error: wrong number of arguments for setfont
) print flush
newpath
10 ux 10 uy moveto
(Hello World)   show stroke 
newpath
10 ux 10 uy moveto
(Hello World)   show stroke 
newpath
10 ux 10 uy moveto
(Hello World)  56 rotate   show stroke  -56 rotate 
(error: text: wrong number of arguments
) print flush
(error: text: wrong number of arguments
) print flush
(error: bad curve definition, wrong number of args
) print flush
newpath
10 ux 310 uy moveto
10 ux 320 uy 110 ux 310 uy 110 ux 320 uy curveto
110 ux 330 uy 210 ux 330 uy 210 ux 320 uy curveto stroke
(error: bad curvextend definition, wrong number of args
) print flush
';
}

