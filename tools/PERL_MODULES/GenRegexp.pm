package GenRegexp;

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


#use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use DataFiles;


sub new {
    my $self  = {};
    
    my $df = DataFiles->new;

    $self->{SEQ1} = undef; 
    $self->{PROGRAM} = $df->get('GENREGEXP');
    $self->{TABLE} = Table->new;
    $self->{TWOSTRAND} = 1;

    bless($self);           # but see below
    return $self;
}


sub setVerbose {

    my ($self, $i) = @_;
    $self->{VERBOSE}    = $i;
    
}


sub setTwoStrands {

    my ($self, $i) = @_;
    $self->{TWOSTRAND}    = $i;
    
}

sub setReldist {

    my ($self, $i) = @_;
    $self->{RELDIST}    = $i;
    
}


sub setSimple{
  my ($self, $i) = @_;
  $self->{SIMPLE} = $i;
}

# 
#
sub run {
    
    my ($self, $re, $fasta) = @_;

    my $t1 = Sets::getTempFile("/tmp/file1");
    
    my $s_todo1 = "$self->{PROGRAM} -re \"$re\" -fastafile $fasta ";
    
    if (defined($self->{RELDIST})) {
	$s_todo1 .= " -reldist $self->{RELDIST} ";
    }
    
    if (defined($self->{TWOSTRAND})) {
	$s_todo1 .= " -twostrand $self->{TWOSTRAND} ";
    }

    if (defined($self->{SIMPLE})) {
	$s_todo1 .= " -simple $self->{SIMPLE} ";
    }

    $s_todo1 .= "> $t1";


    #print "$s_todo1\n";
    print "Scanning genome 1 .. " if ($self->{VERBOSE});
    system($s_todo1);  

    
    $self->{TABLE}->loadFile($t1);
    
    unlink $t1;
}


sub getTable {
    my ($self)=@_;

    return $self->{TABLE};
}


sub getNbMatches {
    my ($self)=@_;

    return $self->{TABLE}->getNbRows();
    
}

sub getMatchingSequences {
  my ($self)=@_;
  
  my @a = ();
  foreach my $r ( @{ $self->{TABLE}->getArray() } ) {
    push @a, $r->[0] if ($r->[1] == 1);
  }
  return \@a;
  
}




1;

