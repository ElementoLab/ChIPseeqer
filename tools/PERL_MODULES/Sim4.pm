package Sim4;

use Sets;
use strict;




sub new {
    my ($self) = {};
   
    $self->{sim4bin}     = "/home/elemento/PERL_MODULES/PROGRAMS/SIM4/sim4.2002-03-03/sim4";
    $self->{XXXX}        = undef;
    $self->{UNAME}       = `uname`; $self->{UNAME} =~ s/\n//g;
    $self->{EXONS}       = [];
    $self->{QLEN}        = undef;
    $self->{DLEN}        = undef;
    bless($self);
    return $self;
}




sub set {
    my ($self, $n)       = @_;
    $self->{XXXX}        = undef;
}

sub getExons {
    my ($self)           = @_;
    return $self->{EXONS};
}

 
#
#  run
#
sub run {
    my ($self, $f1, $f2)           = @_;
    my $tmpfile                    = Sets::getTempFile("/tmp/sim4");
    my $todo = "$self->{sim4bin} $f1 $f2 A=0 W=12 K=11 C=11 > $tmpfile";
    system($todo);

#    system("cat $tmpfile");

    # get the data
    my $a_ref = $self->_get_result($tmpfile);
    

    return $a_ref;
}


sub _get_result {
  my ($self, $f)           = @_;

  open IN, $f;
  my @o = ();
  while (my $l = <IN>) {

    chomp $l;

    if ($l =~ /^seq1/) {
      $l =~ /\,\ (\d+)\ bp/;
      $self->{QLEN} = $1;
      next;
    }

    if ($l =~ /^seq2/) {
      $l =~ /\,\ (\d+)\ bp/;
      $self->{DLEN} = $1;
      next;
    }

    next if ($l eq "");
    next if ($l eq "(complement)");

    #print "$l\n";
    my @a = split /\ +/, $l, -1;
    push @o, \@a;

  }
  close IN;

  my $minc = undef;
  my $maxc = undef;
  my $ming = undef;
  my $maxg = undef;

  my $i = 0;
  
  my $avgsim = undef;
  my $sumlen = 0;

  foreach my $e (@o) {
    my @a = split /\-/, $e->[0];
    
    if ($i == 0) {
      $minc = $a[0];
      my ($i1, $i2) = $e->[1] =~ /\((\d+)\-(\d+)\)/;
      $ming = $i1;
    }
    if ($i == scalar(@o)-1) {
      
      $maxc = $a[1];
      my ($i1, $i2) = $e->[1] =~ /\((\d+)\-(\d+)\)/;
      $maxg = $i2;
    }

    my $l = $a[1] - $a[0] + 1; 
    $sumlen += $l;
    my $id = $e->[2]; $id =~ s/\%$//; $id = int($id);

    $avgsim += $l * $id;

    $i++;
  }


  $avgsim = int(0.5+$avgsim / $sumlen);

  #print "minc=$minc, maxc=$maxc\n";
  #print "ming=$ming, maxg=$maxg\n";

  my %H = ( RAW => \@o, MINC => $minc, MAXC => $maxc, MING => $ming, MAXG => $maxg, AVGSIM => $avgsim );

  return \%H;
}



1;

