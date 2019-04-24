#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;

my $l = <STDIN>;
chomp $l;
my @a = split /\t/, $l, -1;
my @todel = ();
my @toshow = ();
my $i = 0;
my $s = shift @a;
foreach my $r (@a) {
  if (Sets::in_array($r, @ARGV)) {
    push @todel, $i; 
    push @toshow, $r;
  } 
  $i++;
}
print "$s\t" . join("\t", @toshow); print "\n";



@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $n = shift @a;
    my @l = ();
    for (my $i=0; $i<scalar(@a); $i++) {
	$c = $a[$i];
	push @l, $c if (Sets::in_array($i, @todel));
	
    }
    print "$n\t" . join ("\t", @l) . "\n";
}




