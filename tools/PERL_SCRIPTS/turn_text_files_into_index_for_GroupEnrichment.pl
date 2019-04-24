use lib qw(/home/elemento/PERL_MODULES);
use Sets;

my %H = ();
foreach my $f (@ARGV) {
    
    my $s   = Sets::readSet($f);
    my $txt = $f;
    $txt =~ s/\.txt//;

    foreach my $e (@$s) {
	push @{ $H{$e} }, $txt;
    }

}


foreach my $k (keys(%H)) {
    print "$k\t"; print join("\t", @{$H{$k}}); print "\n";
}

