use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $n = scalar(@$a_ref);

my $cnt = $a_ref->[$n-1]->[1];
my $jjj = $a_ref->[$n-1]->[0];

my $slice = int(0.5+$cnt / $ARGV[1]);


my @A = ( [0, 0] );
for (my $i=1; $i<$n; $i++) {
    
    if (int($a_ref->[$i]->[1] / $slice) != int($a_ref->[$i-1]->[1] / $slice)) {
	my @o = ($i, $a_ref->[$i]->[1]); 

	push @A, \@o;
    }
    
} 


$A[$#A]->[0] = $jjj;

$n = scalar(@A);
for (my $i=1; $i<$n; $i++) {
    my $t = $A[$i-1][0]+1; $t = 0 if ($i == 1); 
    print $t . "\t" . $A[$i][0] . "\n";
}
