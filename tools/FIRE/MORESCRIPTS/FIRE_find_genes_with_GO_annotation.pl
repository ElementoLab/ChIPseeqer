use lib "$ENV{FIREDIR}/SCRIPTS";

use Getopt::Long;
use Table;
use Sets;

GetOptions("expfile=s" => \$expfile,
	   "pathway=s" => \$pathway,
	   "bin=s"     => \$bin,
	   "species=s" => \$species);

#
# load pathways
#
my $speciesfile = "$ENV{FIREDIR}/FIRE_DATA/SPECIES_DATA/$species";

#
# load species file
#
my $ta = Table->new;
$ta->loadFile($speciesfile);
my $h_ref = $ta->getIndexKV(0,1);

# load annotation
my $goindexfile = "$ENV{FIREDIR}/$h_ref->{goindex}";
my $gonamesfile = "$ENV{FIREDIR}/$h_ref->{gonames}";
my $descfile    = "$ENV{FIREDIR}/$h_ref->{genedesc}";

# get matching pathway
$ta->loadFile($gonamesfile);
my $a_ref = $ta->getArray();
my %matchinggo = ();
foreach my $r (@$a_ref) {
  if (($r->[1] =~ /$pathway/) || ($r->[0] eq $pathway)) {
    $matchinggo{ $r->[0] } = 1;
  }
}

# get genes in the corresp pathway
$ta->loadFile($goindexfile);
$a_ref = $ta->getArray();

my %genes = ();
foreach my $r (@$a_ref) {
  my $g = shift @$r;
  
  foreach my $c (@$r) {

    if (defined($matchinggo{ $c })) {
      $genes{ $g } = 1; 
    }	

  }
  
}

# load $descfile
$ta->loadFile($descfile);
my $h_ref_desc = $ta->getIndex(0);


# now go thru the expfile 
$ta->loadFile($expfile);
$a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (defined($genes{$r->[0]})) {
    if (!defined($bin) || (defined($bin) && ($r->[1] == $bin))) {
      print "$r->[0]\t$h_ref_desc->{$r->[0]}->[1]\t$h_ref_desc->{$r->[0]}->[2]\n";
    }
  }
}
