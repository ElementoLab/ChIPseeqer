#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my %ADJ = ();
foreach my $r (@$a_ref) {

  next if (
	   ($r->[0] eq "") ||
	   ($r->[0] eq "-") ||
	   ($r->[1] eq "") ||
	   ($r->[1] eq "-") 
	  );

  push @{ $ADJ{ $r->[0] } }, $r->[1] ;

  
  
}


foreach my $g (keys(%ADJ)) {
  print "$g\t" . join("\t", @{ $ADJ{ $g } }) . "\n";
}
