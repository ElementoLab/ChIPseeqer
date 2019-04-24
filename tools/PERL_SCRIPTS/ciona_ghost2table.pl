BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use strict;
my $a_ref = Sets::readSet($ARGV[0]);


my @H            = ();
my %ARRAY        = ();
my @CLUSTERS     = ();
my %NOT_EXAMINED = ();
my %UBIQUITOUS   = ();
my @STAGE_LABELS = ();

foreach my $r (@$a_ref) {

  
  
  my $a_ref_txt = Sets::readSet($r);
  my $n         = scalar(@$a_ref_txt);
  $r =~ s/GHOST\/(\d{5})\.html/$1/;
  #print "$r";

  push @CLUSTERS, $r;

  my $cnt_column = 0;
  for (my $i=0; $i<$n; $i++) {

    if ($a_ref_txt->[ $i ] =~ /\*\*\*\*/) {
      my $st = $a_ref_txt->[ $i   ];  $st =~ s/\n//g;  $st =~ s/<HR>//g;  $st =~ s/\*\*\*\*\* //g; $st =~ s/ \*\*\*\*\*//g;
      $st =~ s/The //g; $st =~ s/ stage//g;

      $STAGE_LABELS[ $cnt_column ] = $st if (!defined($STAGE_LABELS[ $cnt_column ]));

      my $l  = $a_ref_txt->[ $i+1 ];
      $l =~ s/<BR>//g;       $l =~ s/[\r\n]//g; 
      $l =~ s/ +$//g; 
      $l =~ s/^ +//g; 

      #print "<$l>\n";

      if (($l eq "not examined")  || ($l =~ /not clear/))  {

	$NOT_EXAMINED{ $r } [ $cnt_column ] = 1;
      #} else {
      } elsif ($l =~ /ubiquitous/) {
	
	$UBIQUITOUS{ $r } [ $cnt_column ] = 1;
	
      } elsif ($l ne "not detected") {
	
	#print " deal with\n";
	my @a = split /\,/, $l;
	
	foreach my $e (@a) {
	  $H [ $cnt_column ] { $e } = $r;  #, $e if (!Sets::in_array($e, @{ $H [ $cnt_column ] })); 
	  $ARRAY{ $r } [ $cnt_column ] { $e } = 1;
	}

      }

      $cnt_column ++;
      #print "\t$l";
    }
  }
  #print "\n";
}

my @COL_LABELS = ();
my $m_all = scalar(@H);
my $cnt_total_cols = 0;
for (my $i=0; $i<$m_all; $i++) {

  #print "Column $i:\n";

  my @lab = sort(keys(%{ $H[ $i ] }));
  
  #print join("\n ", @lab);
  
  $cnt_total_cols += scalar(@lab);

  $COL_LABELS[ $i ] = \@lab;
  
}


my $n     = scalar(@CLUSTERS);

# traverse the genes ..

print "NONE";
my $cnt_stage = 0;
foreach my $a_ref_type (@COL_LABELS) {
  
  foreach my $t (@$a_ref_type) {
    
    print "\t$STAGE_LABELS[$cnt_stage]/$t";
    
  }
  $cnt_stage ++;


}
print "\n";

      
for (my $i=0; $i<$n; $i++) {
  

  # traverse STATES
  my $cnt_col = 0;
  my $txt     = "";
  my $cnt_X   = 0;
  my $cnt_1   = 0;
  foreach my $a_ref_type (@COL_LABELS) {
    
    # traverse cell types
    foreach my $t (@$a_ref_type) {
      
      #print "$t\n";
      
      #next if ($t eq "not examined");
      
      if (defined($NOT_EXAMINED{ $CLUSTERS[$i] }[ $cnt_col ] )) {
	$txt .= "\t"; $cnt_X ++;
      } elsif (defined($UBIQUITOUS{ $CLUSTERS[$i] }[ $cnt_col ] ) || defined( $ARRAY{ $CLUSTERS[$i] } [ $cnt_col ] { $t } )) {
	$txt .= "\t1"; $cnt_1 ++;
      } else { 
	$txt .= "\t0";
      }
      
    }
    $cnt_col++;
  }
  
  if (($cnt_X != $cnt_total_cols) && ($cnt_1 > 0))  {
    print "$CLUSTERS[$i]$txt\n";
  }

  
  
}






