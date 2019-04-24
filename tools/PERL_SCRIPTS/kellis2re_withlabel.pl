BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Sets;
use Table;

my %H = ( "N" => ".", 
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

           
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $s = $ta->getArray();

foreach my $r (@$s) {

  $r->[1] =~ s/[\r\n]//g; 
  $r->[1] = uc($r->[1]);
  $r->[1] =~ s/^\ +//;
  $r->[1] =~ s/\ +$//;
  my @a = split //, $r->[1], -1;
  print "$r->[0]\t";

  foreach my $t (@a) {
    next if ($t eq ""); 
    die "No match for '$t' in '$r->[1]'\n" if (!defined($H{$t}));
    print "$H{$t}";
  } 
  print "\n";

}
