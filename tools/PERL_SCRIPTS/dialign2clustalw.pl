# usage: perl test_clustalw.pl 992/mavid.mfa > 992/mavid.mfa.pos


use lib qw(/home/elemento/PERL_MODULES);
use Fasta;
use ClustalW;




my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my @a_n = ();
my @a_s = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    push @a_s, $s;
    push @a_n, $n;
    
}

my $w = 200;
my $l = length( $a_s[0] );
my $n = scalar(@a_s);

for (my $i=0; $i<$l-$w; $i+=10) {
  
  my $cnt_ignore = 0;
  foreach my $s (@a_s) {

    my $ss = substr($s, $i, $w);

    #
    # count the number of gaps 
    #
    
    my @a  = split //, $ss;
    my @g  = grep /\-/, @a;
    my $ng = scalar(@g);

    if ($ng > 0.25 * $w) {
      $cnt_ignore ++;
    }


  }
  
  
  if ( ($n - $cnt_ignore) >= 3) {
    
    print "$i\n";
    
  }

  #else {
  #  print substr($a_s[0], $i, $w); print "\n"; 
  #} 
}


#
#  put DroMel first
#

my $i_dmel = undef;
for (my $i=0; $i<$n; $i++) {
  if ($a_n[$i] =~ /dmel/) {
    $i_dmel = $i;
  }
}


my $t = $a_n[ 0 ]; $a_n[ 0 ] = $a_n[ $i_dmel ]; $a_n[ $i_dmel ] = $t;
my $t = $a_s[ 0 ]; $a_s[ 0 ] = $a_s[ $i_dmel ]; $a_s[ $i_dmel ] = $t;


my $cl = ClustalW->new;
$cl->setSequences(\@a_n, \@a_s);
open OUT, ">$ARGV[0].aln";
print OUT $cl->getClustalWformat();
close OUT;
