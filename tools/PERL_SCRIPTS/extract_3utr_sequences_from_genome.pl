#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Table;
use Sequence;
use Fasta;
use Sets;
use Getopt::Long;

use strict;


if (scalar(@ARGV) == 0) {
  print "Usage : extract_3utr_sequences_from_genome.pl --annotation=FILE --genome=FILE --deflen=INT --forcelen=INT --checkmaxlen=INT --minlen=INT --add1=INT\n";
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
my $add1        = 1;
my $lenfile     = undef;
my $noblast     = 0;

GetOptions ('annotation=s'  => \$annotation,
	    'genome=s'      => \$genome,
	    'deflen=s'      => \$deflen,
	    'verbose=s'     => \$verbose,
	    'logfile=s'     => \$logfile,
	    'noblast=s'     => \$noblast,
	    'forcelen=s'    => \$forcelen,
	    'add1=s'        => \$add1,
	    'lenfile=s'     => \$lenfile,
	    'checkmaxlen=s' => \$checkmaxlen,
	    'minlen=s'      => \$minlen);

die "Please format the genome file.\n" if (($noblast == 0) && (! -e "$genome.nhr"));

my $t = Table->new;
$t->loadFile($annotation);

my $s = undef;

if ($noblast == 0) {
  $s = Sequence->new;
  $s->setBlastDB($genome);
}

my $len = $deflen;

my %LEN = ();
my %SEQ  = ();

if (($checkmaxlen == 1) || ($noblast == 1)) {
 
  my $fa = Fasta->new;
  $fa->setFile($genome);

  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $n =~ s/\ .+$//;   
    $LEN{ $n } = length($s);
    if ($noblast == 1) {
      $SEQ{ $n } = $s;
    }
  }

}

my $h_ref_len = undef;
if (defined($lenfile)) {
  $t->loadFile($lenfile);
  $h_ref_len = $t->getIndexKV(0,1);
  
  if (!defined($deflen)) {
    die "--deflen must be defined.\n";
  }

}

if (defined($logfile)) {
  open LOG, ">$logfile" or die "cannot open $logfile.\n";
}

my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );

foreach my $r (@$a_ref) {

  #
  #  case where only the protein sequences are defined
  #
  if (!defined($r->{"START_T"})) {
    $r->{"START_T"} = $r->{"START_P"};
    $r->{"END_T"}   = $r->{"END_P"};    
  }

  my $start = undef;
  my $end   = undef;
  
  #
  #  strand in first 
  #
  if ($r->{STRAND} > 0) {
    
    if (defined($logfile)) {
      my $leny = $r->{END_T} - $r->{END_P};
      print LOG "$r->{ORF}\t$leny\n";
    }
    
    # that's the basic start
    $start = $r->{"END_P"}; if ($add1 == 1) { $start ++ };
    
    # now how much do we add ?

    if (defined($lenfile)) {
      
      

    } elsif (!defined($forcelen)) {
      
      # normal ?
      $end   = $r->{"END_T"};      
      
      # default ?
      if (defined($len) && ($start == $end)) {
	$end += $len;
      }
      
      # but not too much
      if (abs($end - $start) > 5000) {
	$end   = $r->{"END_P"} + 5000;
      }
      
    } else {
      
      # a fixed value
      $end = $start + $forcelen - 1;
      
    }
    
    
  } else {
    
    if (defined($logfile)) {
      my $leny = $r->{START_P} - $r->{START_T};
      print LOG "$r->{ORF}\t$leny\n";
    }
    
    $end   = $r->{"START_P"}; if ($add1 == 1) { $end -- };
    
    if (!defined($forcelen)) {
      
      # normal
      $start = $r->{"START_T"};
      
      # default
      if (defined($len) && ($start == $end)) {
	$start -= $len; 
      }
      
      # but not too much
      if (abs($end - $start) > 5000) {
	$start   = $r->{"START_P"} - 5000;
      }
            
    } else {
      
      # force a fixed value
      $start = $end - $forcelen + 1;
    }
    
    
    # of course, stay on the chr
    if ($start < 0) {
      $start = 1;
    }
    
    
  }
  
  if (abs($start - $end)+1 < $minlen) {
    if ($verbose == 1) {
      my $lll = abs($start - $end);
      print STDERR "$r->{ORF}: Sequence length = $lll (starts at $start, ends at $end), <= $minlen.\n";
    }	
    next;
  }
  
  if (($checkmaxlen == 1) && ($start > $LEN{ $r->{SCAFFOLD} }) && ($end > $LEN{ $r->{SCAFFOLD} })) {
    next;
  }
  
  if (($checkmaxlen == 1) && ($end > $LEN{ $r->{SCAFFOLD} })) {
    $end = $LEN{ $r->{SCAFFOLD} };
  }
  
  # $s->setVerbose(1);

  my $seq = undef;
  if ($noblast == 0) {
    $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
  } else {
    $seq = substr($SEQ{$r->{SCAFFOLD}}, $start - 1, abs($end-$start)+1);
  }
  
  # my $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
  
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
