#!/usr/bin/perl
use strict;
use lib qw(./lib ../lib t/lib);
use Test::Simple tests => 2;
#use Data::Dumper;
use PostScript::Simple;

my $s = new PostScript::Simple(xsize => 200, ysize => 250, eps => 1);

$s->text( 10, 10, 'Hello World' );
$s->text( {align => "right"}, 10, 10, 'Hello World' );
$s->text( 10, 20, '' );
$s->text( 10, 30, "\000" );
$s->text( 10, 40, undef );
$s->text( 10, 50, 'ONE  TWO  THREE~~~~' );
$s->text( {align => "center", rotate => 49}, 40, 80, 'ONE  TWO  THREE~~~~' );
$s->text( 10, 60, join('', map { chr $_ } (0 .. 19)) );
$s->text( 10, 70, join('', map { chr $_ } (20 .. 39)) );
$s->text( 10, 80, join('', map { chr $_ } (120 .. 139)) );
$s->text( 10, 90, join('', map { chr $_ } (140 .. 159)) );
$s->text( 10, 100, '((()))()()()}{}{}][[]]})()})(]' );

ok( length($s->{'pspages'}) eq length(CANNED()) );
ok( $s->{'pspages'} eq CANNED() );

#print STDERR $s->{'pspages'} . "\n";

#print Dumper $s;
#$s->output('text.eps');

sub CANNED {
return 'newpath
10 ux 10 uy moveto
(Hello World)   show stroke 
newpath
10 ux 10 uy moveto
(Hello World)   dup stringwidth pop neg 0 rmoveto show 
newpath
10 ux 20 uy moveto
()   show stroke 
newpath
10 ux 30 uy moveto
(\000)   show stroke 
(error: text: wrong number of arguments
) print flush
newpath
10 ux 50 uy moveto
(ONE  TWO  THREE~~~~)   show stroke 
newpath
40 ux 80 uy moveto
(ONE  TWO  THREE~~~~)  49 rotate   dup stringwidth pop 2 div neg 0 rmoveto show  -49 rotate 
newpath
10 ux 60 uy moveto
(\000\001\002\003\004\005\006\007\010\011\012\013\014\015\016\017\020\021\022\023)   show stroke 
newpath
10 ux 70 uy moveto
(\024\025\026\027\030\031\032\033\034\035\036\037 !"#$%&\')   show stroke 
newpath
10 ux 80 uy moveto
(xyz{|}~\177\200\201\202\203\204\205\206\207\210\211\212\213)   show stroke 
newpath
10 ux 90 uy moveto
(\214\215\216\217\220\221\222\223\224\225\226\227\230\231\232\233\234\235\236\237)   show stroke 
newpath
10 ux 100 uy moveto
(\(\(\(\)\)\)\(\)\(\)\(\)}{}{}][[]]}\)\(\)}\)\(])   show stroke 
';
}
