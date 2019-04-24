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
#
#  build an array $H{ miRNA } = [ kmer1, kmer2 etc ]
#
foreach my $r (@$a_ref_tab) {
    my @b = split /\_/, $r->[2];
    push @{ $H{ $r->[1] } }, @b;
}



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

for (my $i=0; $i<$n; $i++) {
    print scalar(@{$SETS[$i]}) . "\t$KEYS[$i]\n";
}



