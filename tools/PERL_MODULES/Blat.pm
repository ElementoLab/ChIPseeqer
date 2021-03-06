package Blat;
use Sets;
use Sequence;
use strict;
#use DataFiles;

#my $df = DataFiles->new;

sub new {
    my ($self) = {};
    
    $self->{BLAT_BIN}        = "/Users/olivier/bin/blat"; #$df->get("BLASTDIR");
    if (! -e $self->{BLAT_BIN}) { 
      die "Cannot find blat binary\n"; 
    }
    $self->{DATABASE}         = undef;
    $self->{QUERY}            = undef;
    $self->{VERBOSE}          = 0;
    $self->{EVALUE_THRESHOLD} = undef;
    $self->{STORE}            = 1;
    $self->{QUERY_LENGTH}     = undef;
    $self->{HIT_LENGTH}       = undef;
    $self->{HIT_NAME}         = undef;
    $self->{BLAST_PROGRAM}    = "blastn";
    $self->{NB_PROCESSORS}    = undef;
    $self->{CLINE}            = undef;
    $self->{EXITONCRASH}      = 0;
    $self->{CRASHED}          = 0;
    $self->{LOG}              = 0;
    $self->{RAWOUTPUT}        = undef;
    $self->{UNAME}            = `uname`; $self->{UNAME} =~ s/\n//g;
    bless($self);
    return $self;
}

sub DESTROY {
    my ($self) = @_;

    close LOG;
}


sub setExitOnCrash {
  my ($self, $e) = @_;
  $self->{EXITONCRASH} = $e;
}


sub crashed {
  my ($self) = @_;
  return $self->{CRASHED};
}


sub getRawOutput {
  my ($self) = @_;
  return $self->{RAWOUTPUT};
}

sub log {

    my ($self, $f) = @_;
    $self->{LOG} = $f;

    open LOG, ">log.txt" if ($f == 1);
}


sub setNbProcessors {
    my ($self, $f) = @_;
     $self->{NB_PROCESSORS} = $f;
}


sub setBlastProgram {
    
    my ($self, $f) = @_;
    $self->{BLAST_PROGRAM} = $f;
    
}

sub setBlastDir {
    my ($self, $f) = @_;
    $self->{BLAST_DIR} = $f;
}

sub setVerbose {
    my ($self, $f) = @_;
    $self->{VERBOSE} = $f;
}

sub setDatabase {
    my ($self, $f) = @_;
    $self->{DATABASE} = $f;
}


sub setQuery {
    my ($self, $f) = @_;
    $self->{QUERY} = $f;
}

# run blat
sub blat {
    my ($self, $myoutfile) = @_;

    my $q_file = $self->{QUERY};
    my $d_file = $self->{DATABASE};

    my $outfile = undef;
    if (defined($myoutfile)) {
      $outfile = $myoutfile;
    } else {
      $outfile = Sets::getTmpFile("/tmp/blat");
    }

    my $e = "$self->{BLAT_BIN} $d_file $q_file $outfile > /dev/null";

    print LOG "$e\n" if ($self->{LOG} == 1);

    print "$e\n" if ($self->{VERBOSE} == 1);

    system($e) == 0 or die "Cannot execute:\n$e\n";

    $self->{RESFILE} = $outfile;

    #return $self->getExactMatches($outfile);

}

sub cleanup {
  my ($self) = @_;

  unlink $self->{RESFILE};
}

sub getExactMatches {

  my ($self, $myoutfile) = @_;

  my $outfile = undef;

  if (defined($myoutfile)) {
    $outfile = $myoutfile;
  } elsif (defined($self->{RESFILE})) {
    $outfile = $self->{RESFILE};
  } else {
    die "Cannot find any output file for Blat\n";
  }

  open IN, $outfile or die "Cannot open $outfile.\n";

  my $l = <IN>;
  my $l = <IN>;
  my $l = <IN>;
  my $l = <IN>;
  my $l = <IN>;

  my %matches = ();

  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;

    # enforce exact matches
    if (($a[0] == $a[10]) && ($a[1] == 0) && ($a[4] == 0) && ($a[6] == 0)) { 
      
      my $id = $a[9];
      my %M = ( HIT_NAME => $a[13],
		DFROM    => $a[15],
		DTO      => $a[16],
		DSTRAND  => $a[8]  );
      push @{ $matches{$id} }, \%M;
    }
    
  }
  close IN;

  
  return \%matches;
}



1;
