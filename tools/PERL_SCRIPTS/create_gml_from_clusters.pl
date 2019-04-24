use lib qw(/home/olly/PERL_MODULES);
use Table;
use GD;


$im = new GD::Image(1000,1000);

# allocate some colors
$white = $im->colorAllocate(255,255,255);
$black = $im->colorAllocate(0,0,0);
$red = $im->colorAllocate(255,0,0);
$blue = $im->colorAllocate(0,0,255);

# make the background transparent and interlaced
#$im->transparent($white);
#$im->interlaced(true);




# load within-cluster edges
my $ta = Table->new;
my %ADJ = ();
$ta->setDelim(" ");
$ta->setFile($ARGV[0]);
while (my $r = $ta->nextRow()) {
    push @{ $ADJ{ $r->[0] } }, $r->[2];
    push @{ $ADJ{ $r->[2] } }, $r->[0];
} 
$ta->dispose;


# DFS
my @CLUSTERS = ();
my %COL = ();
foreach my $n (keys(%ADJ)) {
    $COL{$n} = 0;
}
my @size = ();
foreach my $n (keys(%ADJ)) {
    if ($COL{$n} == 0) {
        my @a_nodes = ();
        DFS(\%ADJ, $n, \%COL, \@a_nodes);
	#print join("\t", @a_nodes); print "\n";
        #$size [ scalar(@a_nodes) ] ++;
	
	@a_nodes = sort  sortNodesByGroup @a_nodes;
	
	#print join("\t", @a_nodes); <STDIN>;

	push @CLUSTERS, \@a_nodes;

	

    }
} 

my $PI = 3.1415;
my $ID = 1;
my $baseX = 30;
my $baseY = 30;


print "Creator	\"Cytoscape\"
Version	1.0
graph	[
";

#my $nsize = scalar(@size);
for (my $i=0; $i<scalar(@CLUSTERS); $i++) {
    #print "$i\t$size[$i]\n" if (defined($size[$i]));
    #print join("\t", @{ $CLUSTERS[$i] }); print "\n";
    my $m  = scalar( @{ $CLUSTERS[$i] } );

    my $R  = ( 40 * $m ) / ( 2 * $PI );
    
    for (my $j=0; $j<$m; $j++) {
	
	my $x = $baseX + $R + $R * cos( 2 * $PI * ($j / $m));
	my $y = $baseY + $R + $R * sin( 2 * $PI * ($j / $m));



	$INDEX_NODES{ $CLUSTERS[$i]->[$j] } = "-$ID";

	printNode("-$ID", $x, $y, $CLUSTERS[$i]->[$j]);


	$INDEX_NODES_BY_NAME{ "-$ID" } = [ $x, $y, $CLUSTERS[$i]->[$j] ];

	$ID ++;
    
	
    
    }

    $baseX += 2 * $R + 100;
    
    
    #last;
}


$ta->loadFile($ARGV[0]);
my $a_ref_edges_wit = $ta->getArray();
my $MI_wit          = $ta->getColumn(3);

$ta->loadFile($ARGV[1]);
my $a_ref           = $ta->getArray();

my %H = ();
my $a_ref_edges_bet = [];
my $MI_bet = [];
for my $r (@$a_ref) {
    if (!defined($H{$r->[0]}{$r->[2]}) && !defined($H{$r->[2]}{$r->[0]})) {
	push @{ $a_ref_edges_bet }, $r;
	$H{$r->[0]}{$r->[2]} = 1;
	push @{ $MI_bet }, $r->[3];
    }
}

print scalar(@{ $a_ref_edges_bet }); <STDIN>;
print scalar(@{ $a_ref_edges_wit }); <STDIN>;


#my $a_ref_edges_bet = $ta->getArray();
#my $MI_bet          = $ta->getColumn(3);



my @MI_all = ( @$MI_wit, @$MI_bet ); 

my $avg = Sets::average(\@MI_all);
my $std = Sets::stddev (\@MI_all);
my $thr1 = $avg + 1.0 * $std;
my $thr2 = $avg + 1.5 * $std;

print "$thr1\t$thr2\n"; <STDIN>;

my $EID = 1;
foreach my $r (@$a_ref_edges_wit) {
    
    if ($r->[3] >= $thr1) {
	my $lab = "pd";
	if ($r->[3] >= $thr2) {
	    $lab = "pp";
	}

	printEdge("-$EID", $INDEX_NODES{$r->[0]}, $INDEX_NODES{$r->[2]}, $lab);

	
	$im->line($INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[0]} }->[0],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[0]} }->[1],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[2]} }->[0],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[2]} }->[1], ($lab eq "pd"?$blue:$black));
	
    }

    $EID ++;

} 

my %INDEX_EDGES = ();
foreach my $r (@$a_ref_edges_bet) {
    
    next if (defined($INDEX_EDGES{$r->[0]}{$r->[2]}) || defined($INDEX_EDGES{$r->[2]}{$r->[0]}));

    #print join("\t", @$r); print "\n";
    
    #print "$r->[3] >= $thr1 ?\n";

    if ($r->[3] >= $thr1) {

	my $lab = "pd";
	if ($r->[3] >= $thr2) {
	    $lab = "pp";
	}

	printEdge("-$EID", $INDEX_NODES{$r->[0]}, $INDEX_NODES{$r->[2]}, $lab);


	$im->line($INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[0]} }->[0],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[0]} }->[1],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[2]} }->[0],
		  $INDEX_NODES_BY_NAME{ $INDEX_NODES{$r->[2]} }->[1], ($lab eq "pd"?$blue:$black));
	
	#<STDIN>;

	$INDEX_EDGES{$r->[0]}{$r->[2]} = 1;

	$EID ++;
	
    }

    

} #<STDIN>;

	 
print "]\n";

# load between-cluster edges

my @nodes = values(%INDEX_NODES_BY_NAME);

foreach my $n (@nodes) {
    #$im->setAntiAliased($blue);
    #$im->filledArc($x, $y, 30, 30, 0, 360, gdAntiAliased);
    
    $im->filledArc($n->[0], $n->[1], 30, 30, 0, 360, $red);

    #$im->setAntiAliased($black);
    $im->arc($n->[0], $n->[1], 30, 30, 0, 360, $black);	
    $im->string(gdSmallFont, $n->[0] - gdSmallFont->width * length($n->[2]) / 2, $n->[1] - gdSmallFont->height / 2, $n->[2], $black);
    #print "x=$x\ty=$y\n";
}


open OUT, ">toto.png";
binmode OUT;
print OUT $im->png;
close OUT;


sub sortNodesByGroup {
    #my ($a, $b) = @_;

    my ($aa) = $a =~ /\_(\d+)$/;
    my ($bb) = $b =~ /\_(\d+)$/;

    #print "$a $b\n";

    return $aa <=> $bb;

}

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
    my ($edge_id, $edge_tgt, $edge_src, $edge_label) = @_;

print <<EOF
	edge	[
		root_index	$edge_id
		target	$edge_tgt
		source	$edge_src
		label	\"$edge_label\"
	]
EOF

}


sub DFS {
    my ($h_ref_adj, $node, $h_ref_col, $a_ref_nodes) = @_;

    $h_ref_col->{ $node } = 1;
    
    #print "Visiting $node ..\n";

    foreach my $next (@{ $h_ref_adj->{ $node } }) {

        next if (!defined($h_ref_adj->{$next}));


        if ($h_ref_col->{ $next } == 0) {
            DFS($h_ref_adj, $next, $h_ref_col, $a_ref_nodes);
        }
    }
    
    $h_ref_col->{ $node } = 2;

    push @$a_ref_nodes, $node;

}
