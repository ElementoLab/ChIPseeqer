# FASTA iterator
package Fasta;
use Sets;

use strict;
no strict 'refs';

sub new {
    my ($self)  = {};
    $self->{FILE} = undef;
    $self->{HANDLE} = Sets::getRandomString("HANDLE");
    $self->{SEQ}   = undef;
    $self->{NAME}  = undef;
    $self->{END}   = 0;
    bless($self);
    return $self;

}

sub setFile {
    my ($self, $s) = @_;
    $self->{FILE} = $s;
    open $self->{HANDLE}, $self->{FILE} or die "Fasta.pm: Cannot open $self->{FILE}\n";
}

#
#  starts from a text instead of a 
#
sub setText {
    my ($self, $s) = @_;
    $self->{TEXT} = $s;
    
    my @a = split /[\r\n]/, $s;
    
    $self->{TEXT_ARRAY} = \@a;

    $self->{CNT} = 0;
}


sub restart {
    my ($self) = @_;

    $self->dispose;
    open $self->{HANDLE}, $self->{FILE} or die "Cannot open $self->{FILE}\n";
    $self->{END}   = 0;
    

}

sub dispose {
    my ($self) = @_;
    close $self->{HANDLE};
}



sub _readline {
    my ($self) = @_;
    
    if ($self->{FILE}) {
	my $tmp = $self->{HANDLE};
	return <$tmp>;
    } else {
	return $self->{TEXT_ARRAY}->[ $self->{CNT} ++];
    }
    
    
}


# print a Seq
sub writeSeq {
    my ($self, $f, $n, $s) = @_;
    open OUTCH, ">$f" or die "toto\n";
    print OUTCH ">$n\n$s\n";
    close OUTCH;
}

# returns (name, seq)
sub nextSeq {

    
    
    my ($self) = @_;


    #print "getting a new sequence from $self->{FILE}\n";

    if ($self->{END} == 1) {

	return undef;
    }
    

    
    while (1) { 

	my $line = $self->_readline;
	

	#exit;
	# if we are at the end of the file, spit out the sequence, mark finish
	if (!defined($line)) {
	    $self->{END} = 1;
	    
	    if (defined($self->{NAME})) {
		return [ $self->{NAME}, $self->{SEQ} ];
	    } else {
		return undef;
	    }
	}

	chomp $line;
	
	my $trimmed_line = $line;
	$trimmed_line =~ s/\s//g;

	next if (length($trimmed_line) == 0);
	next if ($trimmed_line =~ /^\#/);

	# if the current line starts with >
	if ($line =~ /\>/) {
	    
	    # if the current seq is empty, set the new name
	    if (length($self->{SEQ}) == 0) {
		$self->{NAME} = substr($line, 1);
	    } else {
		# spit out the sequence
		
		my $name = $self->{NAME};
		my $seq  = $self->{SEQ};

		$self->{NAME} = substr($line, 1);
		$self->{SEQ}  = "";
		return [ $name, $seq ];
		
	    }
	} else {
	    #print "Adding $line to current seq\n";
	    $self->{SEQ} .= $line;
	}
    }
    
    
    
}


1;
