use lib qw(/home/olly/PERL_MODULES);

use Sets;

my $files = Sets::getFiles("*.txt");

my %MIRNAS = ();
foreach my $f (@$files) {

    my $s = Sets::readSet($f);

    my $n = Sets::basename($f);
    
    $MIRNAS{ $n } = 1;
    
    foreach my $g (@$s) {
	push @{ $H{$g} }, $n;
	
	
    }
    
    


} 


my @mirnas = keys(%MIRNAS);

print "\t"; print join("\t", @mirnas); print "\n";
foreach my $g (keys(%H)) {
    print "$g";
    foreach my $mi (@mirnas) {
	if (Sets::in_array($mi, @{$H{$g}})) {
	    print "\t1";
	} else {
	    print "\t0";
	}
	    #join("\t", @{$H{$g}}) ."\n";
    }

    print "\n";
}
