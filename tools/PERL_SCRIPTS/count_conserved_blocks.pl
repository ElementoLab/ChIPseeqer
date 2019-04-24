# count occurences of conserved blocks in a set of alignment



use lib qw(/home/olly/PERL_MODULES);




use  ClustalW;
use Sets;


my $a_ref_files = Sets::readSet($ARGV[0]);
my $a_ref_block = Sets::readSet($ARGV[1]);

my %H = ();
foreach my $f (@$a_ref_files) {
    
    open IN, $f;
    my $l = <IN>;
    chomp $l;
    close IN;

    #my $lc = Sets::getComplement($l);

    foreach my $k (@$a_ref_block) {
	my $a_ref_pos = Sets::getREMotifPositions($k, $l);
	$H{$k} += scalar(@$a_ref_pos);
    }


        
    $cnt ++;

    last if ($cnt == 10);
}


foreach my $r (keys(%H)) {
    print "$r\t$H{$r}\n";
}


