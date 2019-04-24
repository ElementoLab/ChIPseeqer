BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Getopt::Long;

my $header = 1;

GetOptions ('table=s'  => \$file,
	    'dict=s'   => \$dict,
	    'col=s'    => \$col,
	    'k=s'      => \$k,
            'header=s' => \$header,
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

if ($header == 1) {
  my $a = shift @$a_ref;
  print join("\t", @$a); print "\n";
}
foreach my $r (@$a_ref) {
    
    if ($h_ref_kv->{ $r->[ $col ] }  ne "") {
        my $m = $h_ref_kv->{ $r->[ $col ] };
        $m =~ s/\ .*$//;
	$r->[ $col ] = $m; #$h_ref_kv->{ $r->[ $col ] };
        print join("\t", @$r); print "\n";
    } else {
        if ($undef eq "inf") {
	  $r->[ $col ] = "inf";
	  print join("\t", @$r); print "\n";

        } elsif ($undef eq 'show') {
	  print join("\t", @$r); print "\n";
	  
  	  #$r->[ $col ] = "inf"; 

        }
    }
    #print join("\t", @$r); print "\n";
}



