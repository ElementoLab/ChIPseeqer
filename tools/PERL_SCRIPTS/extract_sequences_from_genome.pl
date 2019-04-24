#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Table;
use Sequence;
use Sets;
use Fasta;
use Getopt::Long;

use strict;

if (scalar(@ARGV) == 0) {
  print "Usage : extract_upstream_sequences_from_genome.pl --annotation=FILE --genome=FILE --lengthU=INT --lengthD=INT --minlen=INT --checkmaxlen=INT --tssonly=INT --add1=INT\n";
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

GetOptions ('annotation=s'  => \$annotation,
	    'genome=s'      => \$genome,
	    'lengthU=s'     => \$lenU,
	    'lengthD=s'     => \$lenD,
	    'noblast=s'     => \$noblast,
	    'add1=s'        => \$add1,
	    'tssonly=s'     => \$tssonly,
	    'verbose=s'     => \$verbose,
	    'checkmaxlen=s' => \$checkmaxlen,
	    'minlen=s'      => \$minlen);

die "Please format the genome file.\n" if (($noblast == 0) && (! -e "$genome.nhr"));

my $t = Table->new;
$t->loadFile($annotation);


my $s = undef;

if ($noblast == 0) {
  $s = Sequence->new;
  $s->setVerbose(0);
  $s->setBlastDB($genome);
}


my %LEN  = ();
my %SEQ  = ();

if (($checkmaxlen == 1) || ($noblast == 1)) {

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


my $a_ref = $t->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND") );

foreach my $r (@$a_ref) {

    if ($r->{"END_P"} < $r->{"START_P"}) {
	my $tt = $r->{"START_P"};
	$r->{"START_P"} = $r->{"END_P"};
	$r->{"END_P"} = $tt;
    }

    
   
    my $start     = $r->{START_P};
    my $end       = $r->{END_P};

    if ($add1 == 1) { $end -- };
    if ($add1 == 1) { $start ++ };

    next if (($checkmaxlen == 1) && ($start > $LEN{ $r->{SCAFFOLD} }) && ($end > $LEN{ $r->{SCAFFOLD} }));
    
    if (($checkmaxlen == 1) && ($end > $LEN{ $r->{SCAFFOLD} })) {
      $end = $LEN{ $r->{SCAFFOLD} };
    }

    if (defined($minlen) && ( abs($end - $start)+1 < $minlen)) {
      next;
    }
     
    my $seq = undef;
    if ($noblast == 0) {
      $seq = $s->getSequenceFromBlastDB($r->{SCAFFOLD}, $start, $end);
    } else {
      $seq = substr($SEQ{$r->{SCAFFOLD}}, $start - 1, abs($end-$start)+1);
    }
    
    if ($r->{STRAND} < 0) {
	$seq = Sets::getComplement($seq);
    }

    if ($seq && (length($seq) >= $minlen)) {
        print ">$r->{ORF}\n$seq\n\n";
    } else {
      print STDERR "pb with $seq\n";
    }


}
