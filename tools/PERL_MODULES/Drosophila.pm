package Drosophila;


use lib qw(/home/olly/PERL_MODULES);
use Hypergeom;

use Database;

sub new {    
    my $self  = {};
 
    $self->{DB} = Database->new;
    $self->{DB}->connect("FLYDB") or die "Cannot connect\n";
    #$self->{DB}->setVerbose(1);
    $self->{VERBOSE} = 0;
    bless($self);           # but see below
    return $self;    

}

sub getGeneName {

    my ($self, $o) = @_;
    

    $s = "select * from FLY where ID='$o'";
    
    
    my $a = $self->{DB}->queryOneRecord($s);


    return  $a->{NAME};
    
}



sub isOnChromosome {

    my ($self, $o, $c) = @_;
    $s = "select 1 as YES from FLY where ID='$o' and CHROMOSOME like '$c'";
    my $a = $self->{DB}->queryOneRecord($s);
    return  $a->{YES};
    
}

sub isNotOnChromosome {

    my ($self, $o, $c) = @_;
    $s = "select 1 as YES from FLY where ID='$o' and CHROMOSOME not like '$c'";
    my $a = $self->{DB}->queryOneRecord($s);
    return  $a->{YES};
    
}





1;
