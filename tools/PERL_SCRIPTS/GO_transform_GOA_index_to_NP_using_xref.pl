BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;


# load the refLink
$ta->loadFile($ARGV[2]);
my $h_ref_NP_NM = $ta->getIndexKV(3,2);


# load the xref
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my %H = ();

foreach my $r (@$a_ref) {
  
  

  my $n = $r->[1];
  
  next if ($n eq "");

  my $p = $r->[6]; 

  next if ($p eq "");

  #print "'$n'\t'$p'\n";

  my @a = split /\;/, $p;

  foreach my $s (@a) {
    if ($s =~ /\:(.+?)$/) {
      
      my $np = $1;
      
      if (defined($h_ref_NP_NM->{$np})) {
	push @{ $H{ $n } }, $h_ref_NP_NM->{$np}; 
      }
    }
  }

  #print "$n\t" . join("\t", @{ $H{ $n } }); print "\n";
}



#
# finally, go through the index
#

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %G = ();
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  
  foreach my $s (@{ $H{ $n } }) {
    push @{ $G{ $s } }, @$r;
  }

}

foreach my $g (keys(%G)) {
  my $r = Sets::removeDuplicates( $G{ $g } );
  print "$g\t" . join("\t", @$r); print "\n";
}


