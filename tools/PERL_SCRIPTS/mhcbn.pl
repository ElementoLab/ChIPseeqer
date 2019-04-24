use HTML::Form;
use LWP::Simple;
use LWP;

my $HTML = get("http://www.imtech.res.in/raghava/mhcbn/query2.html");


my @forms = HTML::Form->parse($HTML, "http://www.imtech.res.in");

my ($f) = @forms;

@inputs = $f->inputs;

foreach $i (@inputs) {
    print "TYPE=" . $i->type . ", NAME=" .  $i->name . "\n";
    print "POSSIBLE VALUES=\n";
    @vals = $i->possible_values;
    $cnt = 0;
    foreach $v (@vals) {
	print "$cnt\t$v\n";
	$cnt++;
    }
    print "ENTER VALUE [". $f->value($i->name) . "]\n";
    $inp = <STDIN>;
    chomp $inp;
    $f->value($i->name, $vals[$inp]);
}


my $ua = LWP::UserAgent->new;

while (1) {
    $response = $ua->request($f->click);
    $HTML =  $response->content;

    $WEB = $HTML;
    $WEB =~ s/<[^>]*>//g;

    print $WEB;

    my @forms = HTML::Form->parse($HTML, "http://www.imtech.res.in");
    
    print scalar(@forms) . "\n";
    
    last if (scalar(@forms) <= 1);
	
    $f = $forms[1];
    
}
