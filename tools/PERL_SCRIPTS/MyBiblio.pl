use lib qw(/home/olly/PERL_MODULES);
use MyBiblio;


if (scalar(@ARGV) != 2) {

    die "Please provide DB and PMID\n";
}


my $b = MyBiblio->new;
$b->loadDBFile("/home/olly/PERL_MODULES/MyBiblio.db");
$b->setDB($ARGV[0]);

if (!$b->getArticleByPMID($ARGV[1])) {
    die "Article not found ..\n";
    
}



my $id = undef;
if (!($id = $b->db_exists ($ARGV[1])) ) {
    $b->db_add; 
    $b->generateBibTex;

} else {
    print "Article already exists\n";

    # does a link exists ?
    if ($b->db_existsLink($id, $ARGV[0])) {
	print "Link already exists too ..\n";
    } else {
	$b->db_addLink($id, $ARGV[0]);
	print "Link inserted\n";
	$b->generateBibTex;
	
    }	 
	
}

