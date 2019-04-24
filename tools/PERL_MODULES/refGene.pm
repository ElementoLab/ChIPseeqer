package refGene;
use strict;

=head1 NAME

refGene

=head1 SYNOPSIS

  #!/usr/bin/perl
  use lib "$ENV{HOME}/PERL_MODULES";
  use refGene;

  my $re = refGene->new;
  my $a_ref = $re->FindGenesThatOverlapWith("chr17", 7512444, 7531588);
  foreach my $r (@$a_ref) {
    print join("\t", @$r) . "\n";
  }

=cut

use Sets;
use Table;

sub new {
  my ($self) = {};
  $self->{REFGENE} = "$ENV{HOME}/PROGRAMS/ChIPseeqer-1.0/DATA/refGene.txt.25Nov2009";
  $self->{GENES}   = [];
  $self->{EXT}     = 0;

  # load it

  open IN, $self->{REFGENE} or die "Cannot open $self->{REFGENE}\n";

  my %MRNA_LEN = ();
  my %NAME     = ();
  my %CHR      = ();
  my %STRAND   = ();
  my %M_STA    = ();
  my %M_END    = ();
  my %ORF      = ();
  my %E_STA    = ();
  my %E_END    = ();

  while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;
    
    $NAME   { $a[1] } = 1;
    $ORF    { $a[1] } = $a[12];
    $CHR    { $a[1] } = $a[2];
    $STRAND { $a[1] } = ($a[3] eq "+"?1:-1);
    
    # MRNA
    if (defined($M_STA{ $a[1] })) {
      
      # old mRNA length
      my $old_l = $M_END{ $a[1] } - $M_STA{ $a[1] };
      
      # new length
      my $new_l = $a[5]           - $a[4];
      
      if ($new_l > $old_l) {
	# new mRNA larger, update mRNA coordinates
	$M_STA{ $a[1] } = $a[4];
	$M_END{ $a[1] } = $a[5];

	# also update coding sequence coordinates	
	$E_STA{ $a[1] } = $a[6];
	$E_END{ $a[1] } = $a[7];
      }
            
    } else {
      
      # create mRNA coords      
      $M_STA{ $a[1] } = $a[4];
      $M_END{ $a[1] } = $a[5];
      
      # create corresp coding coords
      $E_STA{ $a[1] } = $a[6];
      $E_END{ $a[1] } = $a[7];
    }
    
  }
  
  close IN;
  
  foreach my $g (keys(%NAME)) {
    my @a_tmp = ($g,          # 0
		 $ORF{$g},    # 1
		 $CHR{$g},    # 2
		 $E_STA{$g},  # 3
		 $E_END{$g},  # 4
		 $STRAND{$g}, # 5
		 $M_STA{$g},  # 6
		 $M_END{$g}); # 7

    push @{ $self->{GENES} }, \@a_tmp;
  }
  
  bless($self);
  return $self;
}

sub setExt {
  my ($self, $ext) = @_;
  $self->{EXT} = $ext;
}

sub FindGenesWhosePromotersStartsAt {
  my ($self, $chr, $i) = @_;

  my @a = ();
  foreach my $r (@{$self->{GENES}}) {
    next if ($r->[2] ne $chr);
    my $tss = undef;
    if ($r->[5] == 1) {
      $tss = $r->[6];
    } else {
      $tss = $r->[7];
    }
    if ($tss == $i) { 
      push @a, $r;
    }

  }
  return \@a;
  
}


sub FindGenesThatOverlapWith {
  my ($self, $chr, $i, $j) = @_;
  
  my @a = ();
  foreach my $r (@{$self->{GENES}}) {
    next if ($r->[2] ne $chr);
    if (Sets::sequencesOverlap($i, $j, $r->[6]-$self->{EXT}, $r->[7]+$self->{EXT}) > 0) {
      push @a, $r;
    }

  }
  return \@a;
}


1;
