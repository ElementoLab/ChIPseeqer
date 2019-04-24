package GnuPlot;
use lib qw(/home/olly/PERL_MODULES);

sub new {
    my $self  = {};

    $self->{TYPE} = "png";
    $self->{TITLES} = [];
    $self->{FILES} = [];

    $self->{XLABEL} = "frequency";
    $self->{YLABEL} = "Pearson correlation values";
    
    bless($self);
    return $self;
}

#
# PS or PNG
#
sub setType {
    
    my ($self, $t) = @_;
    
    $self->{TYPE} = $t;
}

#
# title for each curve
#
sub setTitles {
    
    my ($self, $a_ref) = @_;
    
    $self->{TITLES} = $a_ref;
}


#
# title for each curve
#
sub setFiles {
    
    my ($self, $a_ref) = @_;
    
    $self->{FILES} = $a_ref;
}


sub addFile {
    
    my ($self, $file) = @_;
    
    push @{$self->{FILES}}, $file;
    
}

sub addTitle {
    
    my ($self, $t) = @_;
    
    push @{$self->{TITLES}}, $t;
    
}

#
# set the labels 
#
sub setLabels {
    
    my ($self, $lx, $ly) = @_;
    
    $self->{XLABEL} = $lx;
    $self->{YLABEL} = $ly;
        
}


#
# 
#
sub plot {
    
    my ($self, $s_outfile) = @_;
    
    my $s_tmpfile = "/tmp/tmp.plt";
    my $s_plot    = "plot ";
    
    my $i = 0;
    foreach my $f (@{ $self->{FILES} } ) {
	
	#die "$f does not exist ..\n" if (! -e $f);
	my $t = $self->{TITLES}->[$i];

	$s_plot .= "\"$f\" title \"$t\" with lines, ";
	
	$i++;
    }

    chop $s_plot; chop $s_plot;
    $s_plot .= "\n";

    open  OUT, ">$s_tmpfile";
    print OUT "set terminal $self->{TYPE}\n";
    print OUT "set ylabel \"$self->{YLABEL}\"\n";
    print OUT "set xlabel \"$self->{XLABEL}\"\n";
    print OUT "set output \"$s_outfile\"\n";

    print OUT $s_plot;
    close OUT;
    system(("/usr/bin/gnuplot $s_tmpfile"));
    unlink $s_tmpfile;
    
}



1;
