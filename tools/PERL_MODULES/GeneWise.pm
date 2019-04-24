package GeneWise;
use Sets;
use strict;

sub new {
    my ($self) = {};
   
    $ENV{'WISECONFIGDIR'} = "/home/elemento/PERL_MODULES/PROGRAMS/GENEWISE/wise2.2.0/wisecfg/";

    $self->{XXXX}        = undef;
    $self->{UNAME}       = `uname`; $self->{UNAME} =~ s/\n//g;
    $self->{EXONS}       = [];
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
    my $tmpfile                    = Sets::getTempFile("/tmp/genewise");
    my $todo = "/home/elemento/PERL_MODULES/PROGRAMS/GENEWISE/wise2.2.0/src/bin/genewise $f1 $f2 -splice flat -intron tied -genesf 2>/dev/null > $tmpfile";
    system($todo);

#    system("cat $tmpfile");

    # get the data
    my $a_ref = $self->_get_result($tmpfile);

    return $a_ref;
}


sub _get_result {
  my ($self, $f)           = @_;

  open IN, $f;
  my $cnt = 0;
  while ($cnt <= 28) { my $l = <IN>; $cnt++ };
  
  my $do_exit = 0;

  my $s1  = "";
  my $s2  = "";
  my $o_s = undef;
  my $o_e = undef;

  while ($do_exit == 0) {

    $do_exit = 0;
    my @lines = ();
    
    for (my $i=0; $i<8; $i++) {
      my $l = <IN>;
      chomp $l;

      if ($i == 0) {
	my ($s) = $l =~ /^.{20}(.+)$/;
	$s1 .= $s;
      }

      if ($i == 2) {
	my ($s) = $l =~ /^.{20}(.+)$/;
	$s2 .= $s;
      }

      if ($l =~ /^\/\//) {
	$do_exit = 1;
	last;
      }

    }
  }



  
  
  my @a1 = split //, $s1;
  my @a2 = split //, $s2;
  my $n  = scalar(@a1);
  my $ss1 = "";
  my $ss2 = "";
  for (my $i=0; $i<$n; $i++) {
    if (($a1[$i] ne " ") && ($a2[$i] ne " ")) {
      $ss1 .= $a1[$i];
      $ss2 .= $a2[$i];
    }
  }
  
  
  $self->{EXONS} = [];
  
  while (my $l = <IN>) {
    chomp $l;
    if ($l =~ /^Gene (\d+) (\d+)/) {
      $o_s = (!defined($o_s)?$1:$o_s);
      $o_e = ($2>$o_e?$2:$o_e);
    }
    if ($l =~ /  Exon (\d+) (\d+)/) {
      my @a_exon = ($1, $2);
      push @{ $self->{EXONS} }, \@a_exon;
    }
  }
  
  close IN;

  #print "$ss1\n";
  #print "$ss2\n";
  $ss1 =~ s/[\-\ ]//g;
  $ss2 =~ s/[\-\ ]//g;

  return [ $o_s, $o_e, $ss1, $ss2 ];

}



1;

