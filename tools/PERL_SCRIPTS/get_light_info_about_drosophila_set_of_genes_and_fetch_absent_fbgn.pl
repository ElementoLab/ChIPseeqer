#!/usr/bin/perl
#
# this script takes as input a set of FB genes, and look for 1) DB info 2) traductor 
# the assumption is that the traductor
# if nothing is found, then the program tries to get the info from a web site  
#
use lib qw(/home/olly/PERL_MODULES);
use Database;;
use Getopt::Long;
use Sets;
use GO_categories;
use DataFiles;
use Hypergeom;
use Table;
use LWP::Simple;
use LWP;

#
#  load the gene list
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref_genes = $ta->getIndex(0);
my $a_ref = $ta->getColumn(0);

#
#  load the traductor list 
#
$ta->loadFile($ARGV[1]);
my $h_ref_trad = $ta->getIndexKV(0,1);



my $df      = DataFiles->new;

my $db = Database->new();
$db->connect("FLYDB");

my $ntotal = 13522;

#
#   build a set, and get info about it
#
my $set         = Sets::SetToSQLSet($a_ref, "'");
my $sql         = "select * from FLY where ID in $set;";
my $a_ref_genes = $db->queryAllRecordsRef($sql);
my $h_ref       = Sets::SQLRefToSimpleHash($a_ref_genes, "ID", "NAME");

my @server = ("http://fbserver.gen.cam.ac.uk:7081",
	      "http://flybase.bio.indiana.edu",
	      "http://shigen.lab.nig.ac.jp:7081",
	      "http://fly.ebi.ac.uk:7081");

my $txt = "";
foreach my $r (@$a_ref) {

    $r =~ s/gn/GN/g;
   
    #
    #  case 1 : gene is in DB
    #
    if (defined($h_ref->{$r})) {
	print "$r\t$h_ref->{$r}\t\n";
    }

    #
    #  case 2 : trad is available
    #
    elsif (defined($h_ref_trad->{$r})) {
	print "$r\t\t$h_ref_trad->{$r}\n";
    }

    #
    #  case 3 : fetch the FBGN info 
    # 
    else {
		
	my $rr = $r; $rr =~ s/GN/gn/g;
	#print "$rr ?\n";
	my $se = Sets::getRandomElement(\@server);
	my $html = get("$se/.bin/fbidq.html?$rr");
	my ($fbgn) = $html =~ /\<B\>FlyBase\ ID\<\/B\>\<BR\>(FBgn\d+)/;
	
	$fbgn = uc($fbgn);
	print "$r\t\t$fbgn\t$se\n";
	sleep(1+int(rand(5)));

	#my $fbgn = <STDIN>; chomp $fbgn;

	#print "$r\t$h_ref->{$r}\t$fbgn\t$se\n";
	#$txt .= "$r\t$h_ref->{$r}\t$fbgn\t$se\n";
	
    }
    
}




