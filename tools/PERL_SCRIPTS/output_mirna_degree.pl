use lib qw(/home/olly/PERL_MODULES);

use Sets;

my $files = Sets::getFiles("*.txt");

foreach my $f (@$files) {

    my $s = Sets::readSet($f);

    my $n = Sets::basename($f);

    foreach my $g (@$s) {
	push @{ $H{$g} }, $n;
	
	
    }
    
    


} 



foreach my $g (keys(%H)) {
    print "$g\t" . join("\t", @{$H{$g}}) ."\n";
}
