#
#  for each GO term, get its parents
#


use lib qw(/home/elemento/PERL_MODULES);

use GO;
use GO_func;

use strict;
use DataFiles;


my $df = DataFiles->new;

my $go = GO->new;

my $type = $ARGV[0];

if ($type eq "F") {
  $go->loadOntology("./function.ontology");
} elsif ($type eq "P") {
  $go->loadOntology("./process.ontology");
} elsif ($type eq "C") {
  $go->loadOntology("./component.ontology");
} else {
  die "shit..\n";
}


my $a_ref_cat = $go->get_all_categories();

foreach my $r (@$a_ref_cat) {
  my $a_ref_pa = $go->getAllParents($r);
  my @a_allo   = grep !/(^4$|^3673$|^8150$|^3674$|^5575$|^5623$)/, @$a_ref_pa;
  print $r; 
  if (scalar(@a_allo) > 0) {
    print "\t"; print join("\t", @a_allo);
    
  }
  print "\n";
}

