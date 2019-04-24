use lib qw(/home/elemento/PERL_MODULES);

use Sets;
use Table;
# read cluster list
my $h_ref = Sets::getIndex($ARGV[0]);

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();


my %EST = ();
foreach my $r (@$a_ref) {

  my ($id) = $r->[1] =~ /CLSTR(\d{5})r1/;

  if (defined($h_ref->{$id })) {
    push @{ $EST{$id} }, $r->[0];
  }

}

foreach my $k (sort(keys(%EST))) {

  print "$k\t";
  print join("\t", @{ $EST{$k} }); print "\n";


}
