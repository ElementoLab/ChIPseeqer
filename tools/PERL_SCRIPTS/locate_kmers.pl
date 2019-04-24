#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);
use Getopt::Long;
use Table;
use DataFiles;
use strict;

my $utr5file = undef;
my $utr3file = undef;
my $orfsfile = undef;
my $kmerfile = undef;
my $verbose  = 0;
my $bestwins = undef;
GetOptions ('utr5file=s'     =>  \$utr5file,
	    'utr3file=s'     =>  \$utr3file,
	    'orfsfile=s'     =>  \$orfsfile,
	    'kmerfile=s'     =>  \$kmerfile,
	    'bestwins=s'     =>  \$bestwins,
	    'verbose=i'      =>  \$verbose
	    );


if (!$utr5file && !$utr3file && !$orfsfile) {
	die "locate_kmers.pl --utr5file=s --utr3file=s --orfsfile=s --kmerfile=s --bestwins=s --verbose=i\n";
}
	
my $df      = DataFiles->new;
my $s_todo  = undef;


my $ta      = Table->new;
$ta->loadFile($kmerfile);
my $a_ref   = $ta->getArray();

my $h_ref_wins = undef;
if (defined($bestwins)) {
    $ta->loadFile($bestwins);
    $h_ref_wins = $ta->getIndex(0);
}


foreach my $r (@$a_ref) {
    
    $s_todo  = $df->get('KMERS_LOCATE') . " -fasta $utr5file -re $r->[0]";

    if (defined($bestwins)) {
	my $ts = ($h_ref_wins->{$r->[0]}->[2] - $h_ref_wins->{$r->[0]}->[1])/2;
	my $mu = $h_ref_wins->{$r->[0]}->[1] + $ts;
	$s_todo .= " -mu $mu -2s $ts ";
    }
    
    if ($verbose) {
	print "$s_todo\n";
    }
    my ($cnt1, $avg1) = `$s_todo` =~ /^(\d+).+is (\d+)/; 
    $s_todo = $df->get('KMERS_LOCATE')  . " -fasta $orfsfile -re $r->[0]";
    if ($verbose) {
	print "$s_todo\n";
    }
    my ($cnt2, $avg2) = `$s_todo` =~ /^(\d+).+is (\d+)/; 


 #   $s_todo = $df->get('KMERS_LOCATE')  . " -fasta $utr3file -re $r->[0]";
 #   if ($verbose) {
#	print "$s_todo\n";
#    }
    
    #my ($cnt3, $avg3) = `$s_todo` =~ /^(\d+).+is (\d+)/; 
    
    my $cnt3 = "NA";
    
    my $txt   = "nan";


    if ($cnt2 > 0) {
	my $norm = ( $cnt1 * $avg2 ) / ( $cnt2 * $avg1 );
	$txt  = sprintf("%3.2f", $norm);
    }
    

    print "$r->[0]\t$cnt1/$cnt2/$cnt3\t$txt\n";

    
}

