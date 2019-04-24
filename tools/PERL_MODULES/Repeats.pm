package Repeats;

use lib "$ENV{HOME}/PERL_MODULES";

use strict;
use DataFiles;
use File::Basename;
use Sets;


my $df = DataFiles->new;


sub new {
    my ($self)  = {};
    $self->{SEQFILE}        = undef;
    $self->{RESULTS}        = [];
    $self->{TMPFILE}        = 0;
    $self->{PROGRAM}        = "REPEATMASKER";
    $self->{SPECIES}        = undef;
    $self->{TRFMINCOPYNUM}  = undef;
    $self->{TRFMAXUNITSIZE} = undef;
    bless($self);
    return $self;

}

sub setTRFminCopyNumber {
  my ($self, $f) = @_;
  $self->{TRFMINCOPYNUM} = $f;
}

sub setTRFmaxUnitSize {
  my ($self, $f) = @_;
  $self->{TRFMAXUNITSIZE} = $f;
}


sub setSeqFile {
    my ($self, $f) = @_;
    $self->{SEQFILE} = $f;

}


sub setSeq {
    my ($self, $s) = @_;
    
    my $tmpfile = Sets::getTempFile("gogo");
    open TMP, ">$tmpfile";
    print TMP ">TMP\n$s\n";
    close TMP;

    $self->{SEQFILE} = $tmpfile;
    $self->{TMPFILE} = 1;
}

sub setSpecies {
    my ($self, $s) = @_;
    $self->{SPECIES} = $s;
}

sub setProgram {
    my ($self, $s) = @_;
    $self->{PROGRAM} = $s;
}

sub dispose {
    my ($self) = @_;

    if ($self->{TMPFILE} == 1) {
	unlink $self->{SEQFILE};
	unlink "$self->{SEQFILE}.cat";
	unlink "$self->{SEQFILE}.out";
	unlink "$self->{SEQFILE}.tbl";
	unlink "$self->{SEQFILE}.masked";
    }
}

sub process {
    my ($self) = @_;
    
    $self->{RESULTS}   = [];

    
    
    my $todo = undef;

    if ($self->{PROGRAM} eq 'TRF') {

	$todo = $df->get("TRF") . " $self->{SEQFILE} 2 5 5 80 10 30 50 -d > /dev/null";
	#print "$todo\n";

    } elsif ($self->{PROGRAM} eq 'REPEATMASKER') {

	$todo = $df->get("REPEATMASKER") . " $self->{SEQFILE} ";
	if (defined($self->{SPECIES})) {
	    $todo .= " -species $self->{SPECIES}"; 
	}
	$todo .= " &> /dev/null";

    } else {
	die "please specify a program\n";
    }

    system($todo); # == 0 or die "cannot execute TRF\n";

    if ($self->{PROGRAM} eq 'TRF') {
	my $bs   = basename($self->{SEQFILE});
	my $suf  = "$bs.2.5.5.80.10.30.50";
	open IN, "$suf.dat";
	for (my $i=0; $i<15; $i++) {
	    my $l = <IN>;
	}
	
	while (my $l = <IN>) {
	    chomp $l;
	    #print "$l\n";
	    my @a = split / /, $l, -1;
	    
	    # skip if copy numner too low

	    
	    if (defined($self->{TRFMINCOPYNUM}) && ($a[3] < $self->{TRFMINCOPYNUM})) {
	      #print "$a[3] too low, skip\n";
	      next;

	    }

	    # skip if unit size too high
	    if (defined($self->{TRFMAXUNITSIZE}) && ($a[4] > $self->{TRFMAXUNITSIZE})) {
	      #print "$a[4] too high, skip\n";
	      next;
	    }
	    
	    
	    
	    my @a_tmp = ($a[0]-1, $a[1]-1, $a[14]);
	    push @{ $self->{RESULTS} }, \@a_tmp;
	    
	    
	}
	
	close IN;
	
	#print "trying to unlink $suf.dat\n";
	
	unlink "$suf.dat";
	unlink "$suf.1.txt.html";
	unlink "$suf.1.html";

    } elsif ($self->{PROGRAM} eq 'REPEATMASKER') {

	my $outfile   = "$self->{SEQFILE}.cat";
	
	if (! -e $outfile) {
	    die "could not find $outfile\n";
	}
	
	open IN, $outfile;
	while (my $l = <IN>) {
	    last if ($l =~ /sequence/);
	    chomp $l;
	    my @a = split / /, $l, -1;
	    my @a_tmp = ($a[5], $a[6], $a[8]);
	    push @{ $self->{RESULTS} }, \@a_tmp;
	}
	
	#ystem("cat $outfile");
	#system("cat $outfile.out");

	#unlink "$outfile";
	#unlink "$self->{SEQFILE}.out";
	#unlink "$self->{SEQFILE}.tbl";
	#unlink "$self->{SEQFILE}.masked";


    }
}



sub getResults {
    my ($self) = @_;
    return $self->{RESULTS};
}

1;
