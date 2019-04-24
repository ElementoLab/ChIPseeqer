BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use File::Copy;

my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]) or "Cannot load matrix2 file\n";
my $a_ref2 = $ta->getArray();

# get the overlap between the column IDs
my @a1 = @{ $a_ref1->[0] };
my @a2 = @{ $a_ref2->[0] };
shift @a1;
shift @a2;

#my $a_ref_int = Sets::getOverlapSet(\@a1, \@a2);

#print "got " . scalar(@$a_ref_int) . " ov\n";

# matrix 1, what are the columsn to keep ?
my $i = 1;
my %H = ();
foreach my $c (@a1) {
    $H{ $c } = $i;
    $i ++;
}

my @tokeep = ("0");
foreach my $r (@$a_ref_int) {
    push @tokeep, $H{ $r };
}

#print "got " . scalar(@tokeep) . " cols to keeep\n";

copy $ARGV[0],"mm1.txt";


# 

#
# create an index SPECIES => INDEX for the second matrix
#
my $i = 1;
my %H = ();
foreach my $c (@a2) {
    $H{ $c } = $i; #print "pushing $c\n";
    $i ++;
}


open OUT, ">mm2.txt";
# traverse each line


my @new_a1 = ();
foreach my $s (@a1) {
  my $col = $H{ $s };
  if (defined($col)) {
    push @new_a1, $s;
  }
}
print OUT "\t" . join("\t", @new_a1) . "\n";

shift @$a_ref2;
foreach my $r (@$a_ref2) {

    my @row = ();
 
    push @row, $r->[0];
   
    # traverse each line
    foreach my $s (@a1) {
	
	# get the corresponding column
	my $col = $H{ $s };
	
	if (!defined($col)) {
	    #push @row, "inf";
	} else {
	    push @row, $r->[ $col ];
	}
	
    }

    print OUT join("\t", @row); print OUT "\n";
}
close OUT;


