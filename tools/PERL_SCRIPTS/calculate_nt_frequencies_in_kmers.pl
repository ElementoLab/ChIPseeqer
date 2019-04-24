use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getColumn(0);

my $cnt = 0;
foreach my $r (@$a_ref) {
    my @a = split //, $r;
    
    foreach my $s (@a) {
	$H{ $s } ++; $cnt ++;
    }
}


foreach my $k (keys(%H)) {
    print "$k\t" . ($H{$k}/$cnt) . "\n";
}
