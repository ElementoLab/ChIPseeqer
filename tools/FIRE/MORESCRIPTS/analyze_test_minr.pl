use Table;
use Sets;

#my @a = (0.1, 0.5); 

#my @a = (0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,50,100,500,1000,10000,100000);
my @a = (0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,50,100,500,1000,10000,100000);
my $ta = Table->new;

my $a_ref_library = Sets::getFiles("~/DATA/YEASTS/KNOWN_MOTIFS/*.ace");



foreach my $minr (@a) {
  #print "minr = $minr.\n";
  
  my $sum = "TEST/YEAST/GASCH_CLUSTERED/MINR/$minr/Gasch_cinds_OnlyPosCorr.txt.nodups.optim7.summary";

  my %H = ();
  my $cnt_u = 0;

  $ta->loadFile($sum);
  my $a_ref = $ta->getArray();
  my $cnt = 1;
  foreach my $r (@$a_ref) {
    
    my $re = $r->[0];
    
    my $mo = Sets::myre2wm($re);
    open OUT, ">TEST/m$cnt.txt" or die "cannot open mot file for output.\n";
    print OUT "Motif 1\n$mo";
    close OUT;

  
    
    my $best_c = -10000;
    my $best_m = undef;
    foreach my $m (@$a_ref_library) {
      $s_todo = "/home/olly/PERL_MODULES/PROGRAMS/ACE/CompareACE TEST/m$cnt.txt $m";
      $score = `$s_todo`;
      chomp $score;
      if ($score > $best_c) {
	$best_c = $score;
	$best_m = $m;
      }
    }
    
    
    
    
    my $momo = Sets::filename($best_m);
    
    if ($best_c >= 0.8) {
      $H{ $momo } ++;
    } else {
      $cnt_u ++;
    }
    
    #print "$re\t$momo\t$best_c\n";
    
    #system("cat TEST/m$cnt.txt");
    
    $cnt ++;
  }
  
  my $red = scalar(@$a_ref) - ( scalar(keys(%H)) + $cnt_u );

  print "minr = $minr, " . scalar(@$a_ref) . " motifs, " . scalar(keys(%H)) . " unique known motifs, $red redundants, $cnt_u unknown.\n";
  
  

}
