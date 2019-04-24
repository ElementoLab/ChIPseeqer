BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Sets;

my %H = ( "n" => ".", 
	  "A" => "A",
	  "C" => "C",
	  "G" => "G",
	  "T" => "T",
	  "M" => "[AC]",
	  "R" => "[AG]",
	  "W" => "[AT]",
	  "S" => "[CG]",
	  "Y" => "[CT]",
	  "K" => "[GT]",
	  "V" => "[ACG]",
	  "H" => "[ACT]",
	  "B" => "[CGT]",
	  "D" => "[AGT]" );

           

my $s = Sets::readSet($ARGV[0]);

foreach my $r (@$s) {
  
  my @a = split //, $r, -1;

  foreach my $t (@a) {
    print "$H{$t}";
  } 
  print "\n";

}
