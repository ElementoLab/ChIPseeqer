use lib qw(/home/elemento/PERL_MODULES);

use GO;


use strict;
use DataFiles;


#
#  load parents
#
use Table;

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
    # print "$l\n";
    next if ($l =~ /^\!/);

    my @a = split /\t/, $l;

    # GO category
    $a[4] =~ s/GO\://; $a[4] = int($a[4]);

    if (defined($type) && ($type ne $a[8])) {
      next;
    }

    #
    # get the UID
    #
    my $alt      = $a[10];
    my @a_alt    = split /\|/, $alt;
    my @a_alt_cg = grep  /CG/, @a_alt;
    my $uid      = $a_alt_cg[0];
    
    # print "$a[4] on $uid\n";
    
    next if (!defined($uid));

    push @{ $INDEX{ $uid } }, $a[4] if (!Sets::in_array($a[4], @{ $INDEX{ $uid } }));
	    
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




