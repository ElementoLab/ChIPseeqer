#
#  reads in a graph, 
#    # => gene name, adj list of #, randomize gene names, 
#


use lib qw(/home/elemento/PERL_MODULES);

#use Yeast;
use Sets;
use Table;
#use Hypergeom;

my $a_ref = Sets::getArray($ARGV[0]);

my %H = ();  
my $n = 0;
my @ADJ = ();
my @HI = ();

my %INDEX = ();

foreach my $r (@$a_ref) {
    
  next if ( defined($INDEX{ $r->[0] }{ $r->[1] }) || defined($INDEX{ $r->[1] }{ $r->[0] }));
  $INDEX{ $r->[0] }{ $r->[1] } = 1;
  
  if (!defined($H{ $r->[0] })) {
    $H{ $r->[0] } = $n++;
    $HI [ $n-1 ] = $r->[0];
  } 
  
  if (!defined($H{ $r->[1] })) {
    $H{ $r->[1] } = $n++;
    $HI [ $n-1 ] = $r->[1];
  } 
   
  push @{ $ADJ [ $H{ $r->[0] } ] },  $H{ $r->[1] };
  push @{ $ADJ [ $H{ $r->[1] } ] },  $H{ $r->[0] };
   
}


# shuffle names !
#my $a_ref_shu = Sets::shuffle_array( \@HI );



#my @HD = ();
my $cnt = scalar(@ADJ);
print "(mclheader
mcltype matrix
dimensions $cnt" . "x" . "$cnt
)
(mclmatrix
begin
";
my $cnt1 = 0;
foreach my $r1 (@ADJ) {
    
    print "$cnt1\t";
    #print "$a_ref_shu->[$cnt1] => "; 
    foreach my $r2 (@$r1) {
	
	print "$r2 ";
	
	#if (!defined($HD[$cnt1][$r2]) && !defined($HD[$r2][$cnt1])) {
	#    print "$a_ref_shu->[$cnt1]\t$a_ref_shu->[$r2]\n";
	#    $HD[$cnt1][$r2] = 1;
	#} else {
	#print "$a_ref_shu->[$cnt1]\t$a_ref_shu->[$r2]\n";
	#}
    } 
    print "\$\n";
    
    $cnt1++;
}

 print ")\n";
 


open OUT, ">$ARGV[0].dic";
for (my $k=0; $k<scalar(@HI); $k++) {

print OUT "$k\t$HI[$k]\n";


}
close OUT;
