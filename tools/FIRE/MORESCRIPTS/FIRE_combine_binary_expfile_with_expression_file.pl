BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Getopt::Long;
use Table;
use Sets;

my $binaryfile = undef;
my $expfile    = undef;
my $t          = undef;

GetOptions("binaryfile=s"  => \$binaryfile,
	   "expfile=s"     => \$expfile,
	   "t=s"           => \$t);

my $ta = Table->new;	   

# load binaryfile
$ta->loadFile($binaryfile);
my $a_ref = $ta->getArray();

$ta->loadFile($expfile);
my $h_ref = $ta->getIndexKV(0,1);


if (!defined($t)) {
  $t = 0;
} elsif ($t eq "median") {
  my $c = $ta->getColumn(1);
  shift @$c;
  $t = Sets::median($c);
}



print STDERR "Threshold is set to $t\n";



my $r = shift @$a_ref;
print join("\t", @$r) . "\n";
foreach my $r (@$a_ref) {
  
  if (defined($h_ref->{ $r->[0] })) {
    
    if ($r->[1] == 0) {

      print join("\t", @$r) . "\n";

    } else {
      
      if ( $h_ref->{ $r->[0] } > $t) {
	print "$r->[0]\t2\n";
      } else {
	print "$r->[0]\t1\n";
      }
      
    }
    
  }

}

