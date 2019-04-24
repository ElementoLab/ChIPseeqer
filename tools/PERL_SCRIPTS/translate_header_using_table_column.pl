BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Getopt::Long;


#my $undef = "inf";
GetOptions ('table=s'  => \$file,
	    'dict=s'   => \$dict,
	    'col=s'    => \$col,
	    'k=s'      => \$k,
	    'v=s'      => \$v,
	    'undef=s'  => \$undef);



if (!$file) {
    die "Usage : tr.. --table=s --dict=s --col=s --k=s --v=s --undef=s\n";
}


my $ta = Table->new;
$ta->loadFile($file);
my $a_ref = $ta->getArray();

$ta->loadFile($dict);
my $h_ref_kv = $ta->getIndexKV($k, $v);

my $a = shift @$a_ref;

my $n = shift @$a; print "$n";
foreach my $r (@$a) {
  print "\t" . $h_ref_kv->{ $r };
}
print "\n";

foreach my $r (@$a_ref) {
  print join("\t", @$r); print "\n";
}




