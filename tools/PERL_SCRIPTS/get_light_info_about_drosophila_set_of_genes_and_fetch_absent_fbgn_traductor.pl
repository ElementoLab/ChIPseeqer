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
use LWP::Simple;
use LWP;


my $ta = Table->new;
#my $a_ref = Sets::readSet($ARGV[0]);
$ta->loadFile($ARGV[0]);
my $h_ref_genes = $ta->getIndex(0);
my $a_ref = $ta->getColumn(0);
# table !

my $df      = DataFiles->new;

my $db = Database->new();
$db->connect("FLYDB");

my $ntotal = 13522;

#  get gene info
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
   
    if (defined($h_ref->{$r})) {
	#print "$r\t$h_ref->{$r}\t\n";
	print "$r\t$r\n"; #$h_ref->{$r}\t\n";
	#$txt .= "$r\t$h_ref->{$r}\t\n";
    } elsif ($h_ref_genes->{$r}->[2] ne "") {
	print "$r\t" . $h_ref_genes->{$r}->[2] . "\n";
	#print join("\t", @{ $h_ref_genes->{$r} }) . "\n";
	#$txt .= join("\t", @{ $h_ref_genes->{$r} }) . "\n";
    } else {

	
	
	#my $rr = $r; $rr =~ s/GN/gn/g;
	#print "$rr ?\n";
	#my $se = Sets::getRandomElement(\@server);
	#my $html = get("$se/.bin/fbidq.html?$rr");
	#my ($fbgn) = $html =~ /\<B\>FlyBase\ ID\<\/B\>\<BR\>(FBgn\d+)/;
	
	#$fbgn = uc($fbgn);
	
	#my $fbgn = <STDIN>; chomp $fbgn;

	#print "$r\t$h_ref->{$r}\t$fbgn\t$se\n";
	#$txt .= "$r\t$h_ref->{$r}\t$fbgn\t$se\n";
	#sleep(1+int(rand(5)));
    }
    
}




