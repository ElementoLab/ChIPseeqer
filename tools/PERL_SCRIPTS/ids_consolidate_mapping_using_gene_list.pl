BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $a_ref_genes = Sets::getIndex($ARGV[1]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my @a = split /[\ \t]+/, $r->[1];
  foreach my $s (@a) {
    $s =~ s/\.\d+//;	
    if (defined($a_ref_genes->{$s})) {
      print "$r->[0]\t$s\n";
      last;  # one a gene is found, stop ... 
    }
    
  }
}

