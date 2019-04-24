BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
my $l = <IN>;
print $l;
while (my $l = <IN>) {
  chomp $l;
  if ($l eq "") {
    print "$l\n";
  } else {	

    if (!defined($w)) {
      my ($nt) = $l =~ /^(.+?\ +)/;
      $w = length($nt);
      print "$w\n";
      $w  = undef;
    }

  }
}

close IN;
