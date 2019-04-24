use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Log;


my $lo = Log->new;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_int = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $h_ref_dic = $ta->getIndexKV(0, 1);

# read the interactions
my %HH = ();


foreach my $i (@$a_ref_int) {
  if (defined($h_ref_dic->{ $i->[0] }) && defined($h_ref_dic->{ $i->[1] })) {
    
    my $r1 = $h_ref_dic->{ $i->[0] };
    my $r2 = $h_ref_dic->{ $i->[1] };
    $r1 =~ s/\-.+$//;
    $r2 =~ s/\-.+$//;
    
    if (!defined($HH{ $i->[0] }{ $i->[1] }) && !defined($HH{ $i->[1] }{ $i->[0] })) {
      print "$r1\t$r2\t$i->[2]\t$i->[3]\n";
      $HH{ $i->[0] }{ $i->[1] } = 1;
      $HH{ $i->[1] }{ $i->[0] } = 1;
    }

  } else {
    $lo->log("$i->[0]\t$i->[1]\ could not be added to the list because ..\n");
  }

}
# read the id to cg dic
