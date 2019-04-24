use HTML::Form;
use LWP::Simple;
use LWP;

$HTML = get("http://www.yeastgenome.org");

@forms = HTML::Form->parse($HTML, "http://www.yeastgenome.org");

print scalar(@forms) . "\n";

foreach $f (@forms) {
    
    print "FORM\n";
    
    @inputs = $f->inputs;

    foreach $i (@inputs) {
	print $i->type . " " .  $i->name . "\n";
	
    }
    
    <STDIN>;
    
    $f->value("query", "GCN4");
    
    my $ua = LWP::UserAgent->new;

    $response = $ua->request($f->click);
    
    #print $response->content;
    
    $RES = $response->content;
  
    $RES =~ s/<[^>]*>//g;
    
    print $RES;
  
    #while ($RES =~ /\<a href="(.+)"\>(.+)\<\/a\>\<\/td\>/g) {
#	print "A=$2, L=$1\n";
 #   }

    
}
