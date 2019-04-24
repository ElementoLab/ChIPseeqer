use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Sequence;
use Fasta;
use strict;
use Repeats;
#use Data::Dumper;

my $fa = Fasta->new;

my $se = Sequence->new;


my $ta = Table->new;


#
#  load the exons
#

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %CHR            = ();
my %EXONS          = ();
my %BOUNDARIES     = ();
my %BOUNDARIES_CDS = ();
my %STRAND         = ();

foreach my $r (@$a_ref) {
    my @a_exon = ($r->[2], $r->[3]);
    push @{ $EXONS{ $r->[1] } }, \@a_exon;
    $CHR{ $r->[1] }             = $r->[0]; 
    $BOUNDARIES{ $r->[1] }->[0] = (defined($BOUNDARIES{ $r->[1] }->[0])?Sets::min($BOUNDARIES{ $r->[1] }->[0], $r->[2]):$r->[2]);
    $BOUNDARIES{ $r->[1] }->[1] = (defined($BOUNDARIES{ $r->[1] }->[1])?Sets::max($BOUNDARIES{ $r->[1] }->[1], $r->[3]):$r->[3]);
    $STRAND{ $r->[1] }          = $r->[4];

}


my @keys = keys(%EXONS);


my $wmax    = 5000;

my $cnt     = 0;
foreach my $k1 (@keys) {
    
    $cnt ++;

    #next if ($cnt < 1056);
    
    #next if ($k1 ne $ARGV[1]);

    #print "--> $k1 $STRAND{$k1}\n";

    # 
    #  get the coding and surrounding regions
    #

    my $w = ($BOUNDARIES{ $k1 }->[0] - $wmax < 0?$BOUNDARIES{ $k1 }->[0]-1:$wmax);


    my $start = Sets::max(0, $BOUNDARIES{ $k1 }->[0] - $w);
    my $end   = $BOUNDARIES{ $k1 }->[0]; 

    $se->setBlastDB($ARGV[1]);

    #print "getting $CHR{ $k1}, $start, $end\n";
    
    my $seq_masked = $se->getSequenceFromBlastDB($CHR{ $k1}, $start, $end);
    
    next if (!$seq_masked);



    if ($STRAND{$k1} == -1) {
	$seq_masked = getSpeComplement($seq_masked);
    }


    print ">$k1\n$seq_masked\n\n";

    
}


sub getSpeComplement {
    
    my ($str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'T';
	} elsif ($s[$i] eq 'T') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} elsif ($s[$i] eq '-') {
	    $d = '-';
	} else {
	    $d = $s[$i];
	}

	$c .= $d;
    }
    
    return $c;
   
}
