use lib "$ENV{HOME}/PERL_MODULES";
use Table;
use Sets;
use strict;

my $PI = 3.1415;
my $ID = 1;
my $baseX = 30;
my $baseY = 30;


#my $listfile = "/Users/olivier/PROGRAMS/PLASMODIUM/PBM/PWMs/RAW/PBMCORES/AP2_PBM_list";

my $listfile1  = "/Users/olivier/PROGRAMS/PLASMODIUM/REFDATA/sortedAP2s_v2.txt";
my $listfile2 = "/Users/olivier/PROGRAMS/PLASMODIUM/REFDATA/otherAP2s_v2.txt";

my $a_ref1    = Sets::readSet($listfile1);
my $a_ref2    = Sets::readSet($listfile2);
print STDERR "Got " . scalar(@$a_ref1) . " and " . scalar(@$a_ref2) . "\n";
my $a_ref     = [];
@$a_ref       = (@$a_ref1, @$a_ref2);

my %HASEDGE = ();
foreach my $ap2 (@$a_ref) {
  my $myap2 = $ap2;
  if ($ap2 =~ /MAL8P1.153/) {
    $myap2 = "MAL8P1153";
  }
  my $tgtfiles = Sets::getFiles("TARGETS_W/*$myap2*.txt | grep -v orthologs");
  foreach my $f (@$tgtfiles) {
    #print "$ap2\t$f\n";
    open IN, $f;
    while (my $l = <IN>) {
      chomp $l;
      my @a = split /\t/, $l, -1;
      if (Sets::in_array($a[0], @$a_ref)) {
	#print "$ap2\t$a[0]\n";
	$HASEDGE{$ap2}{$a[0]} = 1;
      }
    }
    close IN;


  }	
  #print "\n";
}

#exit;

my $m        = @$a_ref;
my $m1       = @$a_ref1;
my $m2       = @$a_ref2;

print "Creator	\"Cytoscape\"
Version	1.0
graph	[
";


my %IDXNODES = ();
my $R  = ( 100 * $m1 ) / ( 2 * $PI );

# first circle
my $cntnode = 0;
for (my $j=0; $j<$m1; $j++) {
  my $x = $baseX + $R + $R * cos( 2 * $PI * ($j / $m1) - $PI / 2);
  my $y = $baseY + $R + $R * sin( 2 * $PI * ($j / $m1) - $PI / 2);
  printNode("-$ID", $x, $y, $a_ref1->[$j]);
  $IDXNODES{ $cntnode++ } = "-$ID";
  $ID ++;
}

# second circle
for (my $j=0; $j<$m2; $j++) {
  my $x = $baseX + $R + ($R/4) * cos( 2 * $PI * ($j / $m2) - $PI / 2);
  my $y = $baseY + $R + ($R/4) * sin( 2 * $PI * ($j / $m2) - $PI / 2);
  printNode("-$ID", $x, $y, $a_ref2->[$j]);
  $IDXNODES{ $cntnode++ } = "-$ID";
  $ID ++;
}



my $EID = 1;
for (my $j=0; $j<$m; $j++) {
  for (my $k=0; $k<$m; $k++) {
    my $lab = "pp";
    if (defined($HASEDGE{$a_ref->[$j]}{$a_ref->[$k]})) {
      printEdge("-$EID", $IDXNODES{ $j }, $IDXNODES{ $k }, $lab);
      $EID ++;
    }
  }
} 

	 
print "]\n";



sub printNode {
    my ($node_id, $node_x, $node_y, $node_label) = @_;
    
    print <<EOF
	node	[
		root_index	$node_id
		id	$node_id
		graphics	[
			x	$node_x
			y	$node_y
			w	30.0
			h	30.0
			fill	\"#ff9999\"
			type	\"ellipse\"
			outline	\"#000000\"
			outline_width	1.0
		]
		label	\"$node_label\"
	]
EOF
}


sub printEdge {
    my ($edge_id, $edge_src, $edge_tgt, $edge_label) = @_;

print <<EOF
	edge	[
		root_index	$edge_id
		source	$edge_src
		target	$edge_tgt
                graphics        [
                        width   1.0
                        fill    \"#000000\"
                        type    \"line\"
                        Line    [
                        ]
                        source_arrow    0
                        target_arrow    3
                ]

		label	\"$edge_label\"
	]
EOF

}


