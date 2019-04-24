BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;
use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $h_ref = $ta->getIndex(0);

my %CHR  = ();
my %CHR2 = ();
shift @$a_ref;
my %H    = ();

foreach my $r (@$a_ref) {
  
  if (defined($H{ "$r->[1]-$r->[2]-N-$r->[4]" }) || defined($H{ "$r->[1]-N-$r->[3]-$r->[4]" })) {

    #print "WAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    next;
  }

  if ($r->[4] == 1) {
    my %h_tmp1 = ('GENE' => $r->[0], 'POS' => $r->[2], 'TYPE' => 'start', 'STRAND' => $r->[4]);
    my %h_tmp2 = ('GENE' => $r->[0], 'POS' => $r->[3], 'TYPE' => 'end',   'STRAND' => $r->[4]);
    push @{ $CHR{ $r->[1] } }, \%h_tmp1;
    push @{ $CHR{ $r->[1] } }, \%h_tmp2;
  } else {
    my %h_tmp1 = ('GENE' => $r->[0], 'POS' => $r->[3], 'TYPE' => 'start',   'STRAND' => $r->[4]);
    my %h_tmp2 = ('GENE' => $r->[0], 'POS' => $r->[2], 'TYPE' => 'end', 'STRAND' => $r->[4]);
    push @{ $CHR{ $r->[1] } }, \%h_tmp1;
    push @{ $CHR{ $r->[1] } }, \%h_tmp2;
  }
  push @{ $CHR2{ $r->[1] } }, $r;
  $H{ "$r->[1]-$r->[2]-N-$r->[4]" } = 1;
  $H{ "$r->[1]-N-$r->[3]-$r->[4]" } = 1;
}

my $lmax = 1000;

foreach my $chr (keys(%CHR)) {
  
  @{ $CHR{ $chr } } = sort { $a->{POS} <=> $b->{ POS } } @{ $CHR{ $chr } } ;

  for (my $i=0; $i<@{ $CHR{ $chr } }; $i++) {
    #print "\t";
    #print $CHR{$chr}->[$i]->{GENE}; print "\t";
    #print $CHR{$chr}->[$i]->{POS }; print "\t";
    #print $CHR{$chr}->[$i]->{TYPE}; print "\t";
    #print $CHR{$chr}->[$i]->{STRAND}; print "\n";
    
    
    my $p = $CHR{$chr}->[$i]->{POS };
    my $s = $CHR{$chr}->[$i]->{STRAND};
    my $g = $CHR{$chr}->[$i]->{GENE  };

    if ($CHR{$chr}->[$i]->{TYPE} eq 'start') {

      # is it contained within another gene ?
      my $inside = 0;
      foreach my $gg (@{ $CHR2{ $chr } }) {

	next if ($g eq $gg->[0]);

	if (($p > $gg->[2]) && ($p < $gg->[3])) {
	  #print "Inside $gg->[0]\n";
	  $inside = 1; last;
	}
      }
      
      if (!$inside) {
	

	
	my $p0 = undef;

	if ($s == 1) {
	  
	  $p --;
	  
	  if ($i == 0) {
	    $p0 = 1;
	  } else {
	    $p0 = $CHR{$chr}->[$i-1]->{POS };	  
	  }

	  my $p0_1k = $p - ($lmax - 1);

	  if ($CHR{$chr}->[$i-1]->{STRAND} == $s) {
	    # print "private";

	  } else {
	    # print "shared";
	    #$p0 = int(0.5 + $p0 + ($p - $p0) / 2);
	    

	  }
	  
	  $p0 = Sets::max($p0, $p0_1k);
	  
	  print "$g\t$chr\t$p0\t$p\t$s\n";

	} else {

	  $p ++;

	  if ($i+1 == @{ $CHR{ $chr } }) {
	    $p0 = $p;
	  } else {
	    $p0 = $CHR{$chr}->[$i+1]->{POS };
	  }

	  my $p0_1k = $p + ($lmax - 1);

	  if ($CHR{$chr}->[$i+1]->{STRAND} == $s) {
	    #print "private";

	  } else {
	    #print "shared";
	    #$p0 = int( 0.5 + $p + ($p0 - $p) / 2 );

	  }

	  $p0 = Sets::min($p0, $p0_1k);
	  
	  print "$g\t$chr\t$p\t$p0\t$s\n";
	}

	
      }
      
    }

    
  }
}
