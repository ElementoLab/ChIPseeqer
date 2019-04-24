package MotifFile;

use lib qw(/home/olly/PERL_MODULES);
use Motif;

sub new {
    my $self  = {};
    $self->{FILE} = undef; # contains all the sites
    $self->{INDEX} = 0;
    $self->{MOTIFS} = [];
    bless($self);           # but see below
    return $self;
}

sub setFile {
    my ($self, $f) = @_;
    
    $self->{FILE} = $f;

    open IN, $f;
    
    my @lines = <IN>;
    my $file  = join "", @lines;
    my @a_parts = split /Motif \d*\n/, $file;

    my $i = 1;
    foreach my $smot (@a_parts) {
	
	next if ($smot eq "");

	#print "mot=$smot()\n";
	
	open TMP1, ">/tmp/tmpi";
	print TMP1 "Motif $i\n$smot\n";
	close TMP1;

	my $mot = Motif->new;
	$mot->readScanACEMotif("/tmp/tmpi");

	#system("cat /tmp/tmpi");
	

	#print $mot->getSimpleMatrix();

	push @{ $self->{MOTIFS} }, \$mot;
	$i++;
	
	#<STDIN>;
	
    }

}


sub nextMotif {
    my ($self) = @_;


    return $self->{MOTIFS}->[ $self->{INDEX}++ ];
    
}


1;
