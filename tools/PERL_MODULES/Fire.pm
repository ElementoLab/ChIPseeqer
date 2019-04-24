package Fire;

use Table;
use strict;

sub loadFireMotifNames {
  my ($file) = @_;
  my $ta = Table->new;
  $ta->loadFile($file);
  my $h_ref = $ta->getIndexKV(0,1);
  return $h_ref;
}


sub loadFireMotifClusters {
  my ($file) = @_;
  my $ta = Table->new;
  $ta->loadFile($file);
  my $h_ref = $ta->getIndexKV(0,1);
  return $h_ref;
}



sub loadFireGOMotifs {
  my ($file) = @_;
  my $ta = Table->new;
  $ta->loadFile($file);
  my $h_ref = $ta->getIndexKV(0,1);
  return $h_ref;
}


sub loadFireMotifRep {
  my ($file) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($file);
  my $a_ref_mo = $ta->getArray();
  my %STAT         = ();
  
  foreach my $r (@$a_ref_mo) {
    $STAT{$r->[0]}->{F} = [];
    $STAT{$r->[0]}->{N} = [];
  }  

  foreach my $r (@$a_ref_mo) {
    #print "Adding $r->[4] and $r->[5] tp $r->[0] at pos $r->[1]\n";
    $STAT{$r->[0]}->{F}->[$r->[1]] = $r->[4];
    $STAT{$r->[0]}->{N}->[$r->[1]] = $r->[5];
  }
    
  return \%STAT;
}

sub loadFireMotifSummary {
  my ($file) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($file);
  my $a_ref_mo = $ta->getArray();
  my %STAT         = ();
  
  my $i = 0;
  foreach my $r (@$a_ref_mo) {
    
    my %a_tmp = (
		 "CNT"    => $i,
		 "RNA"    => $r->[1],
		 "COPIES" => $r->[2],
		 "MI"     => $r->[3],
		 "RANK"   => $r->[4], 
		 "Z"      => $r->[5], 
		 "R"      => $r->[6], 
		 "S"      => $r->[7],
		 "SEED"   => $r->[8],
		 "DIST"   => $r->[9],
		 "ORIE"   => $r->[10],
		 "CONS"   => $r->[11],
		 "NAME"   => undef);
    
    my @clu = (); for (my $i=12; $i<@$r; $i++) { push @clu, $r->[$i]; };
    $a_tmp{CLU} = \@clu;
    $STAT{ $r->[0] }         = \%a_tmp;   
    $i++;
  }
  
  return \%STAT;
  
}




1;
