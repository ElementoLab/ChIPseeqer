BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;

if (@ARGV == 0) {
  die "Args: regions refgene\n";
}

#
# /home/elemento/PEOPLE/WEIMIN/METHYLATION/PLATFORM/prom_coordinates.txt /home/elemento/PEOPLE/WEIMIN/DESIGN/OLIVIER/BCL6_bound_peak_regions.txt
#

#
# load transcripts
#
$ta->loadFile($ARGV[1]);
my $a_ref_transcripts = $ta->getArray();


# load regions
$ta->loadFile($ARGV[0]);

#open STDOUT, ">toto";

my $a_ref_regions = $ta->getArray();
shift @$a_ref_regions;
print STDOUT "PROBESET\tCHR\tSTART\tEND\tTRANSCRIPTS\tGENES\n";
#my $cntreg = 0;
foreach my $m (@$a_ref_regions) {
  
  
  #$cntreg++;
  
  #next if ($cntreg < 23987);
  #print "$cntreg\t$m->[1]\n";
  
  my $ov = 0;
  my @TS = ();
  my @OV = ();
  foreach my $b (@$a_ref_transcripts) {
    
    
    # skip if not same chromosome
    next if ($m->[1] ne $b->[2]);

    #print "\t$b->[1]\n";a`

    # determine tss, depending on strand
    my $tss = undef;
    if ($b->[3] eq '+') {
      $tss = $b->[4];
    } else {
      $tss = $b->[5];
    }
    
    # d = genomic distance between tss and region 
    my $d = undef;
    if (Sets::sequencesOverlap($tss-10, $tss+10, $m->[2], $m->[3])) {
      $d = 0;
    } else {
      $d = Sets::min( abs($m->[2] - $tss), abs($m->[3] - $tss) );
    }

    push @TS, [$d, $b];

    
    if (Sets::sequencesOverlap($m->[2], $m->[3], $b->[4], $b->[5])) {
      push @OV, $b;
    }

  }
  
  # necessary for chrrandom stupid
  next if (@TS == 0);
  
  # sort by ascending order
  @TS = sort { $a->[0] <=> $b->[0] } @TS;


  

  # show all genes with TSS within 250nt
  my $i = 0;
  my @tsnames = ();
  my @genenames = ();
  while ($TS[$i]->[0] <= 250) {
    my $t = $TS[$i]->[1]->[1];
    my $g = $TS[$i]->[1]->[12];
    push @tsnames,   $t if (!Sets::in_array($t, @tsnames));
    push @genenames, $g if (!Sets::in_array($g, @genenames));
    $i++;
  } 

  # add all overlapping genes
  foreach my $t (@OV) {
    push @tsnames,   $t->[1] if (!Sets::in_array($t->[1], @tsnames));
    push @genenames, $t->[12] if (!Sets::in_array($t->[12], @genenames));
  }
  
  my @ids = @tsnames;

  if (@ids > 0) {
    print STDOUT "$m->[0]\t$m->[1]\t$m->[2]\t$m->[3]\t" . join(",", @ids) . "\t" . join(",", @genenames) . "\n";
  } else {
    
    my $d = $TS[0]->[0]; 
    my $g = $TS[0]->[1]->[12];
    my $t = $TS[0]->[1]->[1];
    my $i = $TS[0]->[1]->[4];
    my $j = $TS[0]->[1]->[5];
    if ($d < 5000) { # rescue TSS within 2000
      print STDOUT "$m->[0]\t$m->[1]\t$m->[2]\t$m->[3]\t$t\n";
    } else {
      #print "$m->[0]\t$m->[1]:$m->[2]-$m->[3] Closest gene is $g/$t, d=$d ($i,$j)\n";
      print STDOUT "$m->[0]\t$m->[1]\t$m->[2]\t$m->[3]\tCLOSEST_$t\_$d" . "bp\n";
    }	
    
    #STDOUT->autoflush(1);
  }
  
}


#close STDOUT;
