BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Table;
use strict;
use Sequence;


my $se = Sequence->new;
$se->setBlastDB($ARGV[1]);

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();



my %IND = ();

foreach my $r (@$a_ref) {

    if ($r->[2] eq "CDS") {
	my ($ci) = $r->[8] =~ /\"(.+?)\"/;

	my $exon = $se->getSequenceFromBlastDB(lc($r->[0]), $r->[3], $r->[4]);

	if ($exon) { 
	  
	  if ($r->[6] eq '-') {
	    $exon = Sets::getComplement($exon);
	    $IND{ $ci } = $exon . $IND{ $ci };
	  } else {
	    $IND{ $ci } .= $exon;
	  }
	  
	}
    } 

}



foreach my $r (keys(%IND)) {
    print ">$r\n$IND{$r}\n\n";
}
