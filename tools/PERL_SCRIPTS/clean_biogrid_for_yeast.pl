BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {

  next if ( ($r->[6] =~ /Phenotyp/) ||
	    ($r->[6] =~ /Synthetic/) ||
	    ($r->[6] =~ /Dosage/) ||	    
	    ($r->[6] =~ /RNA/) );

  if ($r->[9] eq "4932") {
    print join("\t", @$r) . "\n";
  }

}

