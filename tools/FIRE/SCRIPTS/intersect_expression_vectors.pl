use lib "$ENV{FIREDIR}/SCRIPTS";


use Table;
use Sets;
use strict;

if (@ARGV == 0) {
  die "Usage: perl intersect_expression_vectors.pl -expfile1 -expfile2 -outfile\n";
}

my $expfile1         = Sets::get_parameter(\@ARGV, "-expfile1");
my $expfile2         = Sets::get_parameter(\@ARGV, "-expfile2");
my $outfile          = Sets::get_parameter(\@ARGV, "-outfile");

my $ta = Table->new;
$ta->loadFile($expfile1);
my $a_ref = $ta->getArray();


$ta->loadFile($expfile2);
my $h_ref = $ta->getIndex(0);

open OUT, ">$outfile" or die "Cannot open $outfile\n";

foreach my $r (@$a_ref) {
  if (defined($h_ref->{ $r->[0] })) {
    print OUT join("\t", @$r); print OUT "\n";
  }
}

close OUT;


