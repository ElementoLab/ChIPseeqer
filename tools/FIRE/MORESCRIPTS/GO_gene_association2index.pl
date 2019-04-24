BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use GO;


use strict;
use DataFiles;


#
#  load parents
#
use Table;
use Sets;

if (@ARGV == 0) {
  die "Usage: perl GO_gene_association2index.pl simplif_GO_annot GO_parents\n"; 
}

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_parents = $ta->getIndexShifted(0);


open IN, $ARGV[0] or die "Cannot open tab file\n";

my $type = $ARGV[2];


#
#  get all the GO categories associated with a single gene
#
my %INDEX = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  push @{ $INDEX{ $a[0] } }, $a[1] if (!Sets::in_array($a[1], @{ $INDEX{ $a[0] } }));
}


#
#  for each gene, get all the OID associated, and retrieve the parents
#
foreach my $uid (keys(%INDEX)) {
  
  #
  # for each uid, get the parents
  #
  
  my $a_ref_go      = $INDEX{ $uid };

  print "$uid";

  my @o = ();
  foreach my $g (@$a_ref_go) {
    push @o, $g if (!Sets::in_array($g, @o));
    my $a_ref_parents = $h_ref_parents->{ $g };
    foreach my $p (@$a_ref_parents) {
      push @o, $p if (!Sets::in_array($p, @o));
    }
  }
  
  if (scalar(@o) > 0) {
    print "\t" . join("\t", @o); 
  }
  print "\n";
}




