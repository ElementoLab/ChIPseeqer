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

  if (defined($H{ $i->[0] }) && defined($H{ $i->[1] })) {
    
    my $r1 = $H{ $i->[0] };
    my $r2 = $H{ $i->[1] };

    if (!defined($HH{ $i->[0] }{ $i->[1] }{ $i->[3] }) && !defined($HH{ $i->[1] }{ $i->[0] }{ $i->[3] })) {
      
      print "$r1\t$r2\t$i->[2]\t$i->[3]\n";

      $HH{ $i->[0] }{ $i->[1] }{ $i->[3] } = 1;
      $HH{ $i->[1] }{ $i->[0] }{ $i->[3] } = 1;
      
    }
    
  

  } else {
    #$lo->log("$i->[0]\t$i->[1]\ could not be added to the list because ..\n");
  }

}
# read the id to cg dic
