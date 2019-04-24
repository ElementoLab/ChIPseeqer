package Fire;

use Table;

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


sub loadFireMotifSummary {
  my ($file) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($file);
  my $a_ref_mo = $ta->getArray();
  my %STAT         = ();
  
  foreach my $r (@$a_ref_mo) {
    
    my %a_tmp = ( "RNA"    => $r->[1],
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
    
    $STAT{ $r->[0] }         = \%a_tmp;   
  }
  
  return \%STAT;
  
}




1;
