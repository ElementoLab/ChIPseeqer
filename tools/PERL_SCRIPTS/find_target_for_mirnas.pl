use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;
use strict;
use Getopt::Long;

if (scalar(@ARGV) < 2) {
    die "Usage : find_target_for_mirnas.pl --kmers_mirnas --outkmers --setdir --targetdir\n";
}

my $outkmers     = undef;
my $kmers_mirnas = undef;
my $setdir       = undef;  
my $targetdir    = undef;
GetOptions ('kmers_mirnas=s' => \$kmers_mirnas,
	    'outkmers=s'     => \$outkmers,
	    'setdir=s'       => \$setdir,
	    'targetdir=s'    => \$targetdir);



# load the microRNA matching file
my $ta = Table->new;
$ta->loadFile($kmers_mirnas);
my $a_ref_tab = $ta->getArray();

my %H = ();
my @a;
if ($outkmers) {
 
    @a = ();
}
#
#  build an array $H{ miRNA } = [ kmer1, kmer2 etc ]
#
foreach my $r (@$a_ref_tab) {
    if ($r->[3] <= 2) {
	$r->[0] =~ s/\d+\: //;
	$r->[1] =~ s/\d+\://;
	push @{ $H{ $r->[1] } }, $r->[0];
	
	if ($outkmers) {
	    push @a, $r->[0];
	}
	
	#print "pushing $r->[0]\n";
    }
    

}


#print join("\n", keys(%H)); print "\n"; exit();

if ($outkmers) {
    Sets::writeSet(Sets::removeDuplicates(\@a), $outkmers);
}

#my $dir = "ACTUAL_CONSERVED_SETS";
#if (defined($ARGV[1])) {
#	$dir = $ARGV[1];
#}

my @SETS = ();
my @KEYS = ();
foreach my $r (sort(keys(%H))) {

    my @a = @{ $H{ $r } };
    
    # get the union of all k-mers
    my $k = shift @a;
    my $set = Sets::readSet("$setdir/$k.txt");
    foreach $k (@a) {
	my $nextset = Sets::readSet("$setdir/$k.txt");
	$set = Sets::getUnionSet($set, $nextset);
    }
    
    #print "Conserved set for $r has size " . scalar(@$set) . "\n";

    #$Sets::writeSet($set, "TARGETS/$r.txt");

    push @SETS, $set;
    push @KEYS, $r;
}


my $n = scalar(@SETS);

for (my $i=0; $i<$n-1; $i++) {
    
    next if (!defined($SETS[$i]));

    for (my $j=$i+1; $j<$n; $j++) {

	next if (!defined($SETS[$j]));
	#next if (!defined($SETS[$i]));
	
	my $r = Sets::getOverlapSet($SETS[$i], $SETS[$j]);
	
	if ( (scalar(@$r) == scalar(@{$SETS[$i]})) && (scalar(@$r) == scalar(@{$SETS[$j]})) ) {
	    
	    $SETS[$j] = undef;

	    $KEYS[$i] = $KEYS[$i] . "_" . $KEYS[$j];
	
	    #print "collapse $KEYS[$i] and $KEYS[$j]\n";
    
	}

    }

}


for (my $i=0; $i<$n; $i++) {
    if (defined($SETS[$i])) {

	print scalar(@{$SETS[$i]}) . "\t$KEYS[$i]\t";
	
	my @a = split /\_/, $KEYS[$i], -1;
	my @b = ();

	foreach my $r (@a) {
	    push @b, @{$H{$r}};
	}

	
	
	print join("_", @{ Sets::removeDuplicates(\@b) }); print "\n";


	Sets::writeSet($SETS[$i], "$targetdir/$KEYS[$i].txt") if (defined($targetdir));
    }
}



