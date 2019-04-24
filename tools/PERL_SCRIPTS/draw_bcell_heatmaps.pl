#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use strict;

if (@ARGV == 0) {
  die "args: genelist genetosort\n";
}

my %datasets = ( "dalla"        => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DALLAFAVERA/ORF_HSAP/GSE2350_series_matrix-12.txt.clean.rowavg.log.qnorm",
		 "dalla-simple" => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DALLAFAVERA/ORF_HSAP/dalla_favera_simplified.txt",
		 "dalla-normal" => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DALLAFAVERA/ORF_HSAP/dalla_favera_NB_CB_MM.txt",
		 "dfci"         => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DFCI_MONTI/ORF/dfci_dlbcl_Montietal.txt",
		 "staudt"       => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/PNAS_ACGH/staudt_dlbcl_pnas_2008.txt");

my %labels = ( "dfci"    => { 
			     "none"        => "none",
			     "GCB_ABC"     => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DFCI_MONTI/ORF/dfci_dlbcl_Montietal.txt.GCB_ABC_labels",
			     "oxphos_BCR"  => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/DFCI_MONTI/ORF/dfci_dlbcl_Montietal.txt.oxphos_BCR_HR_labels"
			    },
	       "dalla"   => { "none"       => "none" },
	       "dalla-simple"   => { "none"       => "none" },
	       "dalla-normal"   => { "none"       => "none" },
	       "staudt"  => { "none"       => "none", 
			      "GCB_ABC"    => "/Users/olivier/PEOPLE/WEIMIN/DALLAFAVERA/PNAS_ACGH/staudt_dlbcl_pnas_2008.txt.GCB_ABC_labels"  }
	     );

foreach my $k (keys(%datasets)) {

  print "$k ... ";

  # extract
  my $todo = "expression_get_subset_of_expression_profiles.pl $ARGV[0] $datasets{$k} | expression_normalize_rows.pl > $ARGV[0].$k";
  #print "$todo\n";
  system($todo);

  foreach my $lab (keys( %{ $labels{$k} }) ) {

    # draw
    $todo = "draw_expression_heatmap.pl --matrix=$ARGV[0].$k --draw=open --h=30 --xbase=200 ";
    
    if ($k !~ /dalla/) {

      if (($ARGV[1]) && ($lab eq "none")) {
	
	$todo .= " --sortcolsbygene=$ARGV[1] ";

	$todo .= " --sortrowsbycorrel=$ARGV[1] ";

      } elsif ($lab ne "none") {
	
	my $desc = $labels{$k}->{$lab};
	$todo .= " --coldesc=$desc --reordercols=1 --suffix=$lab  --clustrows=1 ";
	
      } elsif ($lab eq "none") {
	
	$todo .= " --clustrows=1 ";

      }
    } else {
      
      $todo .= " --clustrows=1 ";
      
    }
    
    print "$todo\n";
    
    system($todo);
  }
  
  print "Done.\n";
}
