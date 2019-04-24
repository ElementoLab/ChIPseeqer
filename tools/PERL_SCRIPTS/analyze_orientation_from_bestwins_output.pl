# determine which are significant after test
use lib qw(/home/olly/PERL_MODULES);
use Table;


#system("perl -pi -e 's/\ +/\t/g' $ARGV[0]");

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my $n = scalar(@$a_ref);

foreach my $r (@$a_ref) {
    print sprintf("$r->[0]\t$r->[1]\t$r->[2]\t$r->[3]\t$r->[4]\t%3.1e%s\t%3.1e%s\n", 
		  $r->[5], ($r->[5]*$n<=0.05?"**":""),
		  $r->[6], ($r->[6]*$n<=0.05?"**":""))

}
