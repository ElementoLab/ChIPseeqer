BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Table;
use Sequence;
use Fasta;
use Sets;
use Getopt::Long;

use strict;


if (scalar(@ARGV) == 0) {
  print "Usage : extract_5utr_sequences_from_genome.pl --annotation=FILE --genome=FILE --deflen=INT --forcelen=INT --checkmaxlen=INT --minlen=INT --logfile=FILE\n";
  exit(0);
}

my $annotation  = undef;
my $genome      = undef;
my $deflen      = undef;
my $forcelen    = undef;
my $checkmaxlen = 1;
my $minlen      = 5;
my $verbose     = 0;
my $logfile     = undef;

GetOptions ('annotation=s'  => \$annotation,
	    'genome=s'      => \$genome,
	    'deflen=s'      => \$deflen,
	    'verbose=s'     => \$verbose,
	    'logfile=s'     => \$logfile,
	    'forcelen=s'    => \$forcelen,
	    'checkmaxlen=s' => \$checkmaxlen,
	    'minlen=s'      => \$minlen);


my $t = Table->new;
$t->loadFile($annotation);

my $s = Sequence->new;
$s->setBlastDB($genome);
my $len = $deflen;

my %LEN = ();

if ($checkmaxlen == 1) {
 
  my $fa = Fasta->new;
  $fa->setFile($genome);

  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $n =~ s/\ .+$//;   
    $LEN{ $n } = length($s);
  }

}

if (defined($logfile)) {
  open LOG, ">$logfile" or die "cannot open $logfile.\n";
}


my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );

foreach my $r (@$a_ref) {


    my $start = undef;
    my $end   = undef;
	
    if ($verbose == 1) {
      print "$r->{ORF}\t$r->{START_T}\t$r->{START_P}\t$r->{END_P}\t$r->{END_T}\t$r->{STRAND}\n";
    }

    #
    #  strand in first 
    #
    if ($r->{STRAND} > 0) {

      if (defined($logfile)) {
	my $leny = $r->{START_P} - $r->{START_T};
	print LOG "$r->{ORF}\t$leny\n";
      }
      
      # that's the end
      $end = $r->{"START_P"} - 1;

      if ($verbose == 1) {
	print "end is $r->{START_P}\n";
      }

      # now how much do we add ?
      if (!defined($forcelen)) {

	# normal ?
	$start   = $r->{"START_T"};  # should not include the last nt
	if ($verbose == 1) {
	  print "start is $r->{START_T}\n";
	}
	
	# default length ?
	if (defined($deflen) && ($r->{"START_P"} == $r->{"START_T"})) {
	  $start = $end - $deflen;
	}
	
	# but not too much
	if (abs($end - $start) > 5000) {
	  $start   = $end - 5000;
	}
      
      } else {
	
	# a fixed value
	$start = $end - $forcelen;
	
      }
      
      # of course, stay on the chr
      if ($start < 0) {
	$start = 1;
      }
      
    } else {
	
      
      if (defined($logfile)) {
	my $leny = $r->{END_T} - $r->{END_P};
	print LOG "$r->{ORF}\t$leny\n";
      }

      $start   = $r->{"END_P"} + 1;  # this we know

      if (!defined($forcelen)) {
	
	# normal
	$end = $r->{"END_T"};
	
	# default
	if (defined($len) && ($r->{"END_P"} == $r->{"END_T"})) {
	    $end = $start + $len;
	}
	
	# but not too much
	if (abs($end - $start) > 5000) {
	  $end   = $start + 5000;
	}

	
      } else {

	# force a fixed value
	$end = $start + $forcelen;
      }
	
    }

    if (abs($start - $end) < $minlen) {
      if ($verbose == 1) {
	my $lll = abs($start - $end);
	print STDERR "$r->{ORF}: Sequence length = $lll (starts at $start, ends at $end), <= $minlen.\n";
      }	
      next;
    }

    if (($checkmaxlen == 1) && ($start > $LEN{ $r->{SCAFFOLD} }) && ($end > $LEN{ $r->{SCAFFOLD} })) {
      print STDERR "$r->{ORF}: Sequence extends beyond genome.\n";
    }
    
    if (($checkmaxlen == 1) && ($end > $LEN{ $r->{SCAFFOLD} })) {
      $end = $LEN{ $r->{SCAFFOLD} };
    }

    my $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
    
    if ($r->{STRAND} < 0) {
	$seq = Sets::getComplement($seq);
    }

    if ($seq && (length($seq) > 1)) {
        print ">$r->{ORF}\n$seq\n\n";
    }
    
    

}


if (defined($logfile)) {
  close LOG;
}
