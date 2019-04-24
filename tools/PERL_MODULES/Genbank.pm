BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

package Genbank;

use Data::Dumper;

sub new {
  my ($self)  = {};

  $self->{SEQUENCE} = "";
  $self->{FEATURES} = ();
  my %h = ();
  $self->{POINTERS} = \%h;

  bless($self);
  return $self;

}

sub getTaxonomy {
  my ($self) = @_;

  return $self->{FEATURES}->[ $self->{POINTERS}->{TAXONOMY} ]->{ORGANISM};
}


sub getSpecies {
  my ($self) = @_;

  return $self->{FEATURES}->[ $self->{POINTERS}->{TAXONOMY} ]->{SPECIES};
}

sub read {
  my ($self, $f) = @_;

  
  open IN, $f;

  my $lenf         = undef;
  my $current_feat = undef;
  my $current_idx  = undef;
  my $sequence     = undef;
  my $org          = undef;
  my $current_key  = undef;
  my $current_TAG  = undef;

  while (my $l = <IN>) {
    chomp $l;
    
    if ($l =~ /^([A-Z]+)\ +/) {
      $current_TAG = $1;
    }

    if ($current_TAG eq "SOURCE") { 
      
      next if ($l =~ /^SOURCE/);

      if ($l =~ /^  ORGANISM/) {
	
	my $sf2 = substr($l, 12);

	

	my %a = (SPECIES => $sf2, ORGANISM => ""); push @{ $self->{FEATURES} }, \%a;
	$current_idx = scalar( @{ $self->{FEATURES} } ) - 1;
	$self->{POINTERS}->{TAXONOMY} = $current_idx;

      } else {
	
	my $sf2 = substr($l, 12);
	$sf2  = " " . $sf2 if ($self->{FEATURES}->[$current_idx]->{ORGANISM} ne "");
	$self->{FEATURES}->[$current_idx]->{ORGANISM} .= $sf2;
	
      }
      
    } elsif ($l =~ /^(FEATURES\ +?)[A-Z]/) {
      $lenf = length($1);
    
    } elsif ($l =~ /^ORIGIN/) {
      $lenf = undef;

      $sequence = 1;

    } elsif ($current_TAG eq "FEATURES") {

      my $sf1 = substr($l, 0, $lenf);
      my $sf2 = substr($l, $lenf);


      $sf1 =~ s/\ //g;


      if (length($sf1) > 0) {

	$current_feat = $sf1;

	my %a = (NAME => $sf1, POS => $sf2); push @{ $self->{FEATURES} }, \%a;

	$current_idx = scalar( @{ $self->{FEATURES} } ) - 1;
        $current_key = "POS";

      } else {


	if ($sf2 =~ /\/(.+?)\=(.+?)$/) {
      
          my $key = $1;
	  my $val = $2;

	  $val =~ s/^\"//;
	  $val =~ s/\"$//;
          

          $current_key = $key;
          $self->{FEATURES}->[ $current_idx ]->{$key}          = $val;
	  
	} elsif ($sf2 =~ /\/(.+?)$/) { 
	  
        } else {

	  $sf2 =~ s/\"$//;
          $self->{FEATURES}->[ $current_idx ]->{$current_key} .= $sf2;

        }
	
      }  
      
      
    } elsif ($current_TAG eq "ORIGIN") {

      $self->{SEQUENCE} .= substr($l, 10); 


    } elsif ($current_TAG eq "ACCESSION") {
      $self->{ACCESSION} = substr($l, 12);
    }

    

    
  }
  

  close IN;
  

  $self->{SEQUENCE} =~ s/\ //g;

}


sub getAccession {
  my ($self) = @_;
  return $self->{ACCESSION};
}

sub getFeatures {
  my ($self) = @_;

  return $self->{FEATURES};
}


sub getPosFromPOS {
  
  my ($self, $pos) = @_;

  my $beg = undef;
  my $end = undef;

  my $st_global = 1;

  if ($pos =~ /^complement\((.+)\)$/) {
    $st_global = 0;
    $pos = $1;
  }

  if ($pos =~ /^join\((.+)\)$/) {
    $pos = $1;
  }

  my @a = split /\,/, $pos;
  my $seq = "";
  foreach my $pp (@a) {
    
    my $st = 1;
    if ($pp =~ /^complement\((.+)\)$/) {
      $pp = $1;
      $st = 0;
    }

    my ($p1, $p2) = split /\.\./, $pp;

    $p1 =~ s/[\>\<]//g;
    $p2 =~ s/[\>\<]//g;

    #print "$p1-1, $p2-$p1+1\n";

    if ($p1 > $p2) {
      die "ooops\n";
    }

    if (!defined($beg) || ($p1 < $beg)) {
      $beg = $p1;
    }

    if (!defined($end) || ($p2 > $end)) {
      $end = $p2;
    }
    
    
  }

  return [$beg, $end, $st_global];
}

sub getSubseqFromPOS {
  my ($self, $pos) = @_;

  my $st_global = 1;

  if ($pos =~ /^complement\((.+)\)$/) {
    $st_global = 0;
    $pos = $1;
  }

  if ($pos =~ /^join\((.+)\)$/) {
    $pos = $1;
  }

  my @a = split /\,/, $pos;
  my $seq = "";
  foreach my $pp (@a) {
    
    my $st = 1;
    if ($pp =~ /^complement\((.+)\)$/) {
      $pp = $1;
      $st = 0;
    }

    my ($p1, $p2) = split /\.\./, $pp;

    $p1 =~ s/[\>\<]//g;
    $p2 =~ s/[\>\<]//g;
    
    my $subseq = uc( substr( $self->{SEQUENCE}, $p1-1, $p2-$p1+1) ); 

    if ($st == 0) {
      $subseq = Sets::getComplement($subseq);
    }

    $seq .= $subseq;

  }

  $seq = uc($seq);

  if ($st_global == 0) {
    $seq = Sets::getComplement($seq);
  }

  return $seq;

}

sub seq {
  my ($self) = @_;
  
  return $self->{SEQUENCE};
}

1;
