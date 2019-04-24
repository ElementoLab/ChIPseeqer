#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);
use Table;
use Sets;
use Getopt::Long;



my $col      = undef;
my $verbose  = 0;
my $evalfile = undef;
my $cond     = 0;
GetOptions ('col=i'        => \$col,
	    'evalfile=s'   => \$evalfile,
	    'cond=i'       => \$cond,
	    'verbose=i'    => \$verbose);

if (!defined($col) || !defined($evalfile)) {
    die "Usage : eval.. --col= --evalfile=\n";
}

my $t = Table->new;
$t->loadFile($evalfile);
my $a_ref_kmers = $t->getArray;


#$a_ref_kmers = Sets::shuffle_array($a_ref_kmers); 

#$t->loadFile($ARGV[1]);
#my $h_ref_kmers = $t->getIndex(0);
my $n  = scalar(@$a_ref_kmers);

if ($verbose == 1) {
    for (my $i=0; $i<$n; $i++) {
	print "\"" . $a_ref_kmers->[$i]->[$col] . "\"";
	print "\n";
	
    }
    exit;
}

my $i  = 1; #scalar(@$a_ref_kmers);
my $nb = 0;
my $w  = 100;

for (my $i=0; $i<$n-$w; $i++) {
    $nb = 0;
    for (my $j=$i; $j<$i+$w; $j++) {
	if ($cond == 0) {
	    $nb++ if (Sets::trim($a_ref_kmers->[$j]->[$col]) ne "");
	} else {
	    my @a = split /\ /, $a_ref_kmers->[$j]->[$col];
	    $nb++ if ($a[0] > 1); 
	}
    }
    
    
    #print $h_ref_kmers->{$r->[1]}->[4] . "\t$nb\n";
    
    $nb = $nb / $w;

    print "$i\t$nb\n";

    

    #$i++;


}
