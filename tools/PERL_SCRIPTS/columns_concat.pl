#!/usr/bin/perl
@t = ();
while ($l = <STDIN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my @l = ();
  foreach my $c (@ARGV) {
    if ($c =~ /\-/) {
      my @b = split /\-/, $c;
      my @d = ();
      foreach my $e (@b) {
	push @d, $a[$e];
      }
      push @l, join("-", @d);
    } else {
      if ($a[$c] eq '+') {
	$a[$c] = 1;       
      }
      if ($a[$c] eq '-') {
	$a[$c] = -1;       
      }
      push @l, $a[$c];   
    }
  }
  print join ("\t", @l) . "\n";
}




