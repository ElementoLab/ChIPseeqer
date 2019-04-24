# simple sql query 


#!/usr/bin/perl
# takes as input new col, + STDIN ORF tab feature
use lib qw(/home/olly/PERL_MODULES);
use Database;;
use Getopt::Long;
use Sets;
use GO_categories;
use DataFiles;
use Hypergeom;
use Table;
use Getopt::Long;
use Insitus;


my $b_desc = 1;
my $b_func = 1;
my $b_exon = 1;
my $b_orth = 1;
my $b_corr = 0;
my $b_tata = 1;

my $genefile = undef;
my $gotable   = undef;
my $ntotal  = undef;
my $html   = 0;
my $maternal  = undef;
my $b_insitus = 1;

my $p_global  = 1e-5;
my $original = 0;

my $desctable = "FLY";
GetOptions ('desc=s'        => \$b_desc,
	    'func=s'        => \$b_func,
	    'exon=s'        => \$b_exon,
	    'tata=s'        => \$b_tata,
	    'orth=s'        => \$b_orth,
	    'corr=s'        => \$b_corr,
	    'insitus=s'     => \$b_insitus,
	    'gotable=s'       => \$gotable,
	    'desctable=s'    => \$desctable,
	    'ntotal=i'       => \$ntotal,
	    'maternal=i'       => \$maternal,
	    'dupl=i'       => \$b_dupl,
            'original=s'       => \$original,

	    'genefile=s'    => \$genefile);
if (!$genefile) {
	die "Please provide at least a --genefile=set.txt\n";

}
my $ta = Table->new;
$ta->loadFile($genefile);
my $a_ref = $ta->getColumn(0);
my $n  = scalar(@$a_ref);
my $h_ref_genes = $ta->getIndex(0);

#my $a_ref = Sets::readSet($ARGV[0]);


my $df      = DataFiles->new;


my $db = Database->new();
$db->connect("FLYDB");

if (!defined($ntotal)) {
    $ntotal = 13522;
}

my $goc = GO_categories->new;
$goc->setVerbose(0);
$goc->setID($df->get("USER"), $df->get("PASS"), $df->get("HOST"));
$goc->setSpecies("DROSOPHILA");
$goc->setOriginal($original);
$goc->setMaxCategory(1000);
$goc->setMinCategory(5);

if (defined($gotable)) {
    $goc->setAnnotationTable($gotable);
}
 
$goc->setTotalNbORFS($ntotal);
my $h_ref_cat = $goc->getGOCategoriesNames();
$goc->setPvalueThreshold($p_global);
#$goc->setBonferroniCorrection(1);



#if ($b_desc) {

    

#  get gene info
my $set         = Sets::SetToSQLSet($a_ref, "'");
my $sql         = "select * from $desctable where ID in $set;";
my $a_ref_genes = $db->queryAllRecordsRef($sql);
my $h_ref_genes_sql = Sets::SQLRefToIndex($a_ref_genes, "ID");
my @H = (); 

my $ORTH  = 0;
my $HUMAN_ORTH  = 0;
my $YEAST_ORTH  = 0;

my $MATER = 0;
my $DUPL  = 0;
my $TATA  = 0;

printHTML("<table>");
foreach my $g (@$a_ref) {
    my $a_ref_cat = $goc->getGOCategories($g);
    
    my $r = $h_ref_genes_sql->{$g};
    
    $H[ $r->{NBEXONS} ] ++ if defined($r->{NBEXONS});
    
    if ($b_desc) {

	printHTML("<tr>\n");

	printHTML("<td valign=top><tt>"); print $g;                             printHTML("</td>"); print "\t";
	printHTML("<td valign=top><tt>"); print $h_ref_genes->{$g}->[1];        printHTML("</td>"); print "\t";
	printHTML("<td valign=top><tt>"); print $r->{NAME};                     printHTML("</td>"); print "\t";
	printHTML("<td valign=top><tt>"); print "EX=$r->{NBEXONS}";             printHTML("</td>"); print "\t";
	printHTML("<td valign=top><tt>"); print "EH=$r->{ESSENTIAL_HEMOCYTES}"; printHTML("</td>"); print "\t";
	printHTML("<td valign=top><tt>"); print "CE=$r->{ELEGANS_ORTHOLOG}";    printHTML("</td>"); print "\t";  
	printHTML("<td valign=top><tt>"); print "MA=$r->{MATERNAL}";            printHTML("</td>"); print "\t";  
	printHTML("<td valign=top><tt>"); print "DU=$r->{DUPLICATE}";           printHTML("</td>"); print "\n";  

	printHTML("<td valign=top><tt>");
	foreach my $c (@$a_ref_cat) {
	    print "\t"; 
	    print "$g\t" . $h_ref_cat->{$c}->{NAME}; 
	    print "\n";
	    printHTML("<br>");
	}
	printHTML("</td>");

	printHTML("</tr>");
	
	printHTML("<tr><td colspan=7>&nbsp;</td></tr>\n");

	
	print "\n";
    }

    if ($r->{ELEGANS_ORTHOLOG} == 1) {
	$ORTH++;
    }

    if ($r->{HUMAN_ORTHOLOGS} == 1) {
	$HUMAN_ORTH++;
    }

    if ($r->{YEAST_ORTHOLOGS} == 1) {
	$YEAST_ORTH++;
    }

    if ($r->{MATERNAL} == 1) {
	$MATER++;
    }
    
    if ($r->{DUPLICATE} == 1) {
	$DUPL++;
    }

    if ($r->{TATA} == 1) {
	$TATA++;
    }
}
printHTML("</table>\n");


my $nn = scalar(@$a_ref);

print "Statistical enrichment for the $nn genes in this set :\n";

if ($b_func) {
    
    print "\n--> Functional enrichment (GO):\n";


    $goc->setORFset($a_ref);
    $goc->setOriginal($original);

    print "1. Biological processes:\n";
    $goc->setAnnotationType("P");
    my $a_ref_func = $goc->getFunctionalContent;
    if (scalar(@$a_ref_func) > 0) {

	foreach my $c (@$a_ref_func) {
	    my $s_func = lc($c->{TEXT});
	    my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	    print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected)\n";
	}

    } else {
	print "Nothing significant\n";
	
    }


    
    print "\n2. Molecular function:\n";
    $goc->setAnnotationType("F");
    my $a_ref_func = $goc->getFunctionalContent;
    if (scalar(@$a_ref_func) > 0) {

	foreach my $c (@$a_ref_func) {
	    my $s_func = lc($c->{TEXT});
	    my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	    print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected)\n";
	}

    } else {
	print "Nothing significant\n";
	
    }

    print "\n3. Cellular localization:\n";
    $goc->setAnnotationType("C");
    my $a_ref_func = $goc->getFunctionalContent;
    if (scalar(@$a_ref_func) > 0) {

	foreach my $c (@$a_ref_func) {
	    my $s_func = lc($c->{TEXT});
	    my $s_prob = sprintf("p=%3.2e", $c->{PVALUE});
	    print "$s_prob\t$s_func\t($c->{OVERLAP} in this set, $c->{EXP} expected)\n";
	}

    } else {
	print "Nothing significant\n";
	
    }


}


if ($b_exon) {
    #  get exon info
    my $sql         = "select NBEXONS, count(*) as COUNT from $desctable group by NBEXONS;";
    my $a_ref_exons = $db->queryAllRecordsRef($sql);
    my $h_ref_exons = Sets::SQLRefToSimpleHash($a_ref_exons, "NBEXONS", "COUNT");
    

    print "\n--> Number of exons / gene:\n";
    
    my $o  = scalar(@H);  
    my $count_exons = 0;
    for (my $i=0; $i<$o; $i++) {
	
	next if (!$H[$i]); 
	my $p1 = Hypergeom::cumbino($H[$i], $n, $h_ref_exons->{$i}/$ntotal);	
	my $p2 = 1.0 - Hypergeom::cumbino($H[$i]+1, $n, $h_ref_exons->{$i}/$ntotal);
	
	#if ( ($p1 < 0.05 / $o) || ($p2 < 0.05 / $o)) {
	if ( ($p1 < $p_global) || ($p2 < $p_global)) {
	    $p = sprintf("%3.2e/%3.2e", $p1, $p2);
	    my $exp = int(0.5 * $n *  $h_ref_exons->{$i}/$ntotal);
	    print "p=$p\tnb exons = $i\t($H[$i] in this set, $exp expected)\n";
	    $count_exons ++;
	}
    }

    if ($count_exons == 0) {
	print "Nothing significant\n";
    }
    
    #print "\n";
}

if ($b_orth) {
    #
    #  get orth info
    #
    my $sql         = "select count(*) as COUNT from $desctable where ELEGANS_ORTHOLOG = 1;";
    my $h_ref_orth = $db->queryOneRecord($sql);
    
    

    print "\n--> Elegans orthologs ($h_ref_orth->{COUNT} genome-wide):\n";
    
    
    
    my $p1 = Hypergeom::cumbino($ORTH, $n, $h_ref_orth->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($ORTH+1, $n, $h_ref_orth->{COUNT}/$ntotal);
        
    #if (($p1 <= 0.05) || ($p2 < 0.05)) {
    if (($p1 <= $p_global) || ($p2 < $p_global)) {
    
	$p = sprintf("%3.2e/%3.2e", $p1, $p2);
	my $exp = int(0.5 + $n * $h_ref_orth->{COUNT} / $ntotal);
	print "p=$p\t($ORTH in this set, $exp expected)\n";
    } else {
	print "Nothing significant\n";
    }

    
    #
    #  get orth info
    #
    my $sql         = "select count(*) as COUNT from $desctable where HUMAN_ORTHOLOGS = 1;";
    my $h_ref_orth = $db->queryOneRecord($sql);
    
    

    print "\n--> Human orthologs ($h_ref_orth->{COUNT} genome-wide):\n";
    
    
    
    my $p1 = Hypergeom::cumbino($HUMAN_ORTH, $n, $h_ref_orth->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($HUMAN_ORTH+1, $n, $h_ref_orth->{COUNT}/$ntotal);
        
    #if (($p1 <= 0.05) || ($p2 < 0.05)) {
    if (($p1 <= $p_global) || ($p2 < $p_global)) {
	    
	$p = sprintf("%3.2e/%3.2e", $p1, $p2);
	my $exp = int(0.5 + $n * $h_ref_orth->{COUNT} / $ntotal);

	print "p=$p\t($HUMAN_ORTH in this set, $exp expected)\n";
    } else {
	print "Nothing significant\n";
    }


    #
    #  get orth info
    #
    my $sql         = "select count(*) as COUNT from $desctable where YEAST_ORTHOLOGS = 1;";
    my $h_ref_orth = $db->queryOneRecord($sql);
    
    

    print "\n--> Yeast orthologs ($h_ref_orth->{COUNT} genome-wide):\n";
    
    
    
    my $p1 = Hypergeom::cumbino($YEAST_ORTH, $n, $h_ref_orth->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($YEAST_ORTH+1, $n, $h_ref_orth->{COUNT}/$ntotal);
        
    #if (($p1 <= 0.05) || ($p2 < 0.05)) {
    if (($p1 <= $p_global) || ($p2 < $p_global)) {
		
	$p = sprintf("%3.2e/%3.2e", $p1, $p2);
	my $exp = int(0.5 + $n * $h_ref_orth->{COUNT} / $ntotal);

	print "p=$p\t($YEAST_ORTH in this set, $exp expected)\n";
    } else {
	print "Nothing significant\n";
    }


    

}


if ($b_insitus) {

    my $in = Insitus->new;

    $in->setORFset($a_ref);
    $in->setBonferroniCorrection(1);
    $in->setPvalueThreshold($p_global);

    $in->calculateEnrichments();
    
    print "\n--> Developmental stages\n";
    my $a_ref_stages = $in->getStageEnrichments();
    if (scalar(@$a_ref_stages)) {
	foreach my $r1 (@$a_ref_stages) {
	    print sprintf("p=%3.2e\t%s\t(%d in this set, %d expected)\n", $r1->{PVALUE}, $r1->{STAGE}, 
			  $r1->{OVERLAP}, $r1->{EXP});
	}
    } else {
	print "Nothing significant\n";
    }
    
    #print "\n--> Tissues\n";
    #my $a_ref_tissues = $in->getTissueEnrichments();
    #if (scalar(@$a_ref_tissues)) {
#	foreach my $r1 (@$a_ref_tissues) {
#	    print sprintf("p=%3.2e\t%s\t(%d in this set, %d genome-wide)\n", $r1->{PVALUE}, $r1->{TISSUE}, 
#			  $r1->{OVERLAP}, $r1->{S1}, $r1->{S2});
#	}
#    } else {
#	print "Nothing significant\n";
#    }
    
    
    print "\n--> Developmental stages / tissues\n";
    my $a_ref_stages_tissues = $in->getStageTissueEnrichments();
    if (scalar(@$a_ref_stages_tissues)) {
	
	my %STAGES = ();
	foreach my $r1 (@$a_ref_stages_tissues) {
	    my ($thisstage) = $r1->{STAGE_TISSUE} =~ /(stage\d+\-\d+)/;
	    my $tmp =  sprintf("p=%3.2e\t%s\t(%d in this set, %d expected)\n", $r1->{PVALUE}, $r1->{STAGE_TISSUE}, 
			  $r1->{OVERLAP}, $r1->{EXP});
	    push @{ $STAGES { $thisstage } }, $tmp;
	}

	my @key_stages = ("stage1-3", "stage4-6", "stage7-8", "stage9-10", "stage11-12", "stage13-16");
	foreach my $st (@key_stages) {
	    #print sprintf("p=%3.2e\t%s\t(%d in this set, %d expected)\n", $r1->{PVALUE}, $r1->{STAGE_TISSUE}, 
	    #		  $r1->{OVERLAP}, $r1->{EXP});
	    next if (scalar(@{ $STAGES { $st } }) == 0);
	    print join("", @{ $STAGES { $st } }); print "\n";
	}

    } else {
	print "Nothing significant\n";
    }
    #print "\n";

}



if ($maternal) {
    #  get orth info
    my $sql         = "select count(*) as COUNT from $desctable where MATERNAL = 1;";
    my $h_ref_mat = $db->queryOneRecord($sql);
    print "\nMaternal genes:\n";
    my $p1 = Hypergeom::cumbino($MATER, $n, $h_ref_mat->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($MATER+1, $n, $h_ref_nat->{COUNT}/$ntotal);
    $p = sprintf("%3.2e/%3.2e", $p1, $p2);
    print "$MATER\tp=$p (genome prop = $h_ref_mat->{COUNT}/$ntotal)\n";

    

}


if ($b_dupl) {
    #  get orth info
    my $sql         = "select count(*) as COUNT from $desctable where DUPLICATE = 1;";
    my $h_ref_dupl = $db->queryOneRecord($sql);
    print "\n--> Duplicate genes ($h_ref_dupl->{COUNT} genome-wide):\n";
    my $p1 = Hypergeom::cumbino($DUPL, $n, $h_ref_dupl->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($DUPL+1, $n, $h_ref_dupl->{COUNT}/$ntotal);
    if (($p1 <= $p_global) || ($p2 < $p_global)) {

	#if (($p1 <= 0.05) || ($p2 < 0.05)) {
	$p = sprintf("%3.2e/%3.2e", $p1, $p2);
	my $exp = int(0.5  + $n * $h_ref_dupl->{COUNT}/$ntotal);
	print "p=$p\t($DUPL in this set, $exp expected)\n";
    } else {
	print "Nothing significant\n";
    }
    

}


if ($b_tata) {
    #  get orth info
    my $sql         = "select count(*) as COUNT from $desctable where TATA = 1;";
    my $h_ref_tata = $db->queryOneRecord($sql);
    print "\n--> consensus TATA-box (TATAAA) containing genes ($h_ref_tata->{COUNT} genome-wide):\n";
    my $p1 = Hypergeom::cumbino($TATA, $n, $h_ref_tata->{COUNT}/$ntotal);
    my $p2 = 1.0-Hypergeom::cumbino($TATA+1, $n, $h_ref_tata->{COUNT}/$ntotal);
    if (($p1 <= $p_global) || ($p2 < $p_global)) {

	#if (($p1 <= 0.05) || ($p2 < 0.05)) {
	$p = sprintf("%3.2e/%3.2e", $p1, $p2);
	my $exp = int(0.5  + $n * $h_ref_tata->{COUNT}/$ntotal);
	print "p=$p\t($TATA in this set, $exp expected)\n";
    } else {
	print "Nothing significant\n";
    }
    

}



if ($b_corr) {
    
    # get the expression profiles


    # calculate the correlation values


    
    
}



sub printHTML {
    my ($x) = @_;

    if ($html) {
	print $x;
    }

}
