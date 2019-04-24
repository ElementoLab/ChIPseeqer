package LynxGet;

sub new {
    my ($self) = {};
    
    
    bless($self);
    return $self;
}

#
#  
#
sub get {
    
    my ($self, $u) = @_;
    
    my $s = `lynx --dump $u`;

    #print $s;

    my @a = split /\n/, $s;

    my @a_links = ();
   
    my $ref = 0;
    foreach my $l (@a) {
	if ($l =~ /\[(\d+)\](.+)$/) {
	    $a_links[$1][0] = $2;
	}
	
	$ref = 1 if ($l =~ /References/);

	next if ($ref == 0);
	

	if ($l =~ /\s+(\d+)\.\s(.+)$/) {
	    $a_links[$1][1] = $2;
	}
    }
    return \@a_links;
}

sub getDump {

     my ($self, $u, $p) = @_;
    
     my $s = `lynx --dump $u $p`;
     
     
}


sub saveDump {

    my ($self, $file, $u, $p) = @_;

    my $s = $self->getDump($u, $p);

    open OUT, ">$file" or die "Could not open $file ..\n"; 
    print OUT $s;
    close OUT;
}


1;

