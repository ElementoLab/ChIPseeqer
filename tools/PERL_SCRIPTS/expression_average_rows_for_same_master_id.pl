use lib qw(/home/elemento/PERL_MODULES);
use strict;
use Sets;

# read in the MID => probes ids file
use Table;

my %H = ();

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $mid = shift @$r;
  foreach my $s (@$r) {
    push @{ $H{ $s } }, $mid; 
  }
}


# read in expression data
my %ROWS = ();

open IN, $ARGV[1];
my $l = <IN>;
print $l;
my $cnt = 0;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $n = shift @a;
  
  #
  # add the line to all MID the ps belongs to
  #
  
  foreach my $mid (@{ $H{ $n } }) {
    push @{ $ROWS{ $mid } }, \@a;
  }
}
close IN;


my @LL = sort(keys(%ROWS));

foreach my $k (@LL) {
    
    my $n = scalar(@{ $ROWS{ $k } });
    my $l = scalar(@{ $ROWS{ $k }->[0] });
    
    my @exp = ();
    for (my $i=0; $i<$l; $i++) {
	my $cnt = 0;
	my $sum = undef;
	for (my $j=0; $j<$n; $j++) {
	    if ($ROWS{ $k }->[$j]->[$i] ne "") {
		$sum += $ROWS{ $k }->[$j]->[$i];
		$cnt += 1;
	    }
	}
	if ($cnt > 0) {
	    $exp[$i] = int( 0.5 + $sum / $cnt ); #sprintf("%3.1f", $sum / $cnt);
	} else {
	    $exp[$i] = "";
	}

	
    }
    
    print "$k\t"; print join("\t", @exp); print "\n";

}
