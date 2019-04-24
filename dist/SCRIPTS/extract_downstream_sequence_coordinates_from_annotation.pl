#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Table;
use Sets;
use Fasta;
use Getopt::Long;

use strict;

if (scalar(@ARGV) == 0) {
	print "Usage : extract_downstream_sequences_from_genome.pl --annotation=FILE --genome=FILE --lengthU=INT --lengthD=INT --minlen=INT --checkmaxlen=INT --tssonly=INT --shortgenes=INT --add1=INT\nFormat is ORF, SCAFFOLD, START_P, END_P, STRAND, START_T, END_T\n";
	exit(0);
}

my $annotation  = undef;
my $genome      = undef;
my $lenU        = undef;
my $lenD        = undef;
my $checkmaxlen = 1;
my $minlen      = 5;
my $verbose     = 0;
my $tssonly     = 0;
my $add1        = 1;
my $noblast     = 0;
my $shortgenes  = 0;

GetOptions ('annotation=s'  => \$annotation,
'genome=s'      => \$genome,
'lengthU=s'     => \$lenU,
'lengthD=s'     => \$lenD,
'noblast=s'     => \$noblast,
'add1=s'        => \$add1,
'tssonly=s'     => \$tssonly,
'verbose=s'     => \$verbose,
'checkmaxlen=s' => \$checkmaxlen,
'shortgenes=i'	=> \$shortgenes,
'minlen=s'      => \$minlen);



my $t = Table->new;
$t->loadFile($annotation);


my $s = undef;


my %LEN  = ();
my %SEQ  = ();

if ($checkmaxlen == 1) { 
	
	my $fa = Fasta->new;
	$fa->setFile($genome);
	
	while (my $a_ref = $fa->nextSeq()) {
		my ($n, $s) = @$a_ref;
		$n =~ s/\ .+$//;
		$LEN{$n} = length($s);
		if ($noblast == 1) {
			$SEQ{ $n } = $s;
		}
	}
	
}


my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );

foreach my $r (@$a_ref) {
	
    if ($r->{"END_T"} < $r->{"START_T"}) {
		my $tt = $r->{"START_T"};
		$r->{"START_T"} = $r->{"END_T"};
		$r->{"END_T"} = $tt;
    }
	
    if ($r->{"END_P"} < $r->{"START_P"}) {
		my $tt = $r->{"START_P"};
		$r->{"START_P"} = $r->{"END_P"};
		$r->{"END_P"} = $tt;
    }
	
    
    #
    #  case where only the protein sequences are defined
    #
    if (!defined($r->{"START_T"})) {
		$r->{"START_T"} = $r->{"START_P"};
		$r->{"END_T"}   = $r->{"END_P"};
		
    }
    
	
    my $start     = undef;
    my $end       = undef;
	
    
    if ($r->{STRAND} > 0) {
		
		if (($tssonly == 1) && ($r->{START_T} == $r->{START_P})) {
			next;
		}
		
		if($r->{END_T} - $lenU > $r->{START_T}) {
			$start = Sets::max(0, $r->{END_T} - $lenU);
			$end   = $r->{END_T} + $lenD; if ($add1 == 1) { $end -- };
		}
		else {
			if($shortgenes == 1) {
				$start	= $r->{START_T};
				$end	= $r->{END_T} + $lenD;
			}
			else {
				next;
			}
		}
		
    } else {
		
		if (($tssonly == 1) && ($r->{END_P} == $r->{END_T})) {
			next;
		}
		
		if($r->{START_T} + $lenU < $r->{END_T}) {
			$start = Sets::max(0, $r->{START_T} - $lenD); if ($add1 == 1) { $start ++ };
			$end   = $r->{START_T} + $lenU; 
		}
		else {
			if($shortgenes == 1) {
				$start = Sets::max(0, $r->{START_T} - $lenD);
				$end   = $r->{END_T}; 
			}
			else {
				next;
			}
		}
    }
    
	
    next if (($checkmaxlen == 1) && ($start > $LEN{ $r->{SCAFFOLD} }) && ($end > $LEN{ $r->{SCAFFOLD} }));
    
    if (($checkmaxlen == 1) && ($end > $LEN{ $r->{SCAFFOLD} })) {
		$end = $LEN{ $r->{SCAFFOLD} };
    }
	
    if ( abs($end - $start)+1 < $minlen) {
		next;
    }
    
    print "$r->{ORF}\t$r->{SCAFFOLD}\t$start\t$end\t$r->{STRAND}\n";
}
