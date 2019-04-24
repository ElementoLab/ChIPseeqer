package TRF;

use lib "$ENV{HOME}/PERL_MODULES";

use strict;
use DataFiles;
use File::Basename;
use Sets;

sub new {
    my ($self)  = {};
    $self->{SEQFILE}   = undef;
    $self->{RESULTS}   = [];
    $self->{TMPFILE}   = 0;
    bless($self);
    return $self;

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

sub dispose {
    my ($self) = @_;

    if ($self->{TMPFILE} == 1) {
	unlink $self->{SEQFILE};
    }
}

sub process {
    my ($self) = @_;
    
    $self->{RESULTS}   = [];
    
    my $todo = "$ENV{HOME}/PROGRAMS/SOFT/trf404.mac-leopard.exe $self->{SEQFILE} 2 5 5 80 10 30 3 -d > /dev/null";
    
    system($todo); # == 0 or die "cannot execute TRF\n";

    my $bs   = basename($self->{SEQFILE});

    my $suf  = "$bs.2.5.5.80.10.30.3";

    open IN, "$suf.dat";
    for (my $i=0; $i<15; $i++) {
	my $l = <IN>;
    }

    while (my $l = <IN>) {

	#print $l;
	chomp $l;

	my @a = split / /, $l, -1;

	my @a_tmp = ($a[0], $a[1], $a[14]);

	push @{ $self->{RESULTS} }, \@a_tmp;
	
    }

    close IN;

    #print "trying to unlink $suf.dat\n";

    unlink "$suf.dat";
    unlink "$suf.1.txt.html";
    unlink "$suf.1.html";

}



sub getResults {
    my ($self) = @_;
    return $self->{RESULTS};
}

1;
