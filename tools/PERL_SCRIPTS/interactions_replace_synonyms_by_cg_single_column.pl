use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Log;
use strict;

my $lo = Log->new;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_int = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $a_ref_dic = $ta->getArray;

my %H = ();
foreach my $r (@$a_ref_dic) {
  my $n = shift @$r;
  foreach my $s (@$r) {
    $H{ $s } = $n;
  }
  
}

# read the interactions
my %HH = ();
foreach my $i (@$a_ref_int) {

  if (defined($H{ $i->[0] })) { 
    
    my $r1 = $H{ $i->[0] };
      
    print "$r1";
    shift @$i; if (scalar(@$i) > 0) { print "\t" . join("\t", @$i); }
    print "\n";
    

  } else {
    #$lo->log("$i->[0]\t$i->[1]\ could not be added to the list because ..\n");
  }

}
# read the id to cg dic
