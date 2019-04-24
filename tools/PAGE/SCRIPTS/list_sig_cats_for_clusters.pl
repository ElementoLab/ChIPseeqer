use Table;
use Sets;
use Getopt::Long;
use strict;

my $expfile           = undef;
my $max_p             = 0.05 ;

GetOptions ('expfile=s'           => \$expfile,
	    'max_p=s'             => \$max_p) ;

#
# creating the summary file
#
my $file = "$expfile\_PAGE/pvmatrix.txt" ;
my $ta = Table->new;
$ta->loadFile($file);

my $a_ref_M      = $ta->getArray();
# header
my $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
if (!defined($max_p)) {
  $max_p = 0.05 / @$a_ref_H;
}

my %CL ;
my $out = "$expfile\_PAGE/list.txt" ;
open O, "> $out" ;
for (my $i=0; $i<@$a_ref_M; $i++) {
 my $r  = $a_ref_M->[$i];
 my $go = shift @$r;
 $go =~ s/^(.+?)\ //;
 $go = "$go, $1" if ($go ne $1) ;
 
 print "$go\n";
 for (my $j=0; $j<@$r; $j++) {
     my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/;
     if ($lpo<log($max_p)){
	 push(@{$CL{$a_ref_H->[$j]}->{go}}, $go) ;
	 push(@{$CL{$a_ref_H->[$j]}->{pv}}, $lpo) ;
     }
 }
}

for (my $i=0 ; $i<@$a_ref_H ; $i++){
    if (! defined $CL{$a_ref_H->[$i]}){
	print O $a_ref_H->[$i], "\n" ;
	next ;
    }
    my @go = @{$CL{$a_ref_H->[$i]}->{go}} ;
    my @pv = @{$CL{$a_ref_H->[$i]}->{pv}} ;

    my @temp = () ;
    for (my $j=0 ; $j<@go ; $j++){
	$temp[$j]->{go} = $go[$j] ;
	$temp[$j]->{pv} = $pv[$j] ;
    }

    @temp = sort {$a->{pv} <=> $b->{pv}} (@temp) ;
    print O $a_ref_H->[$i] ;
    for (my $j=0 ; $j<@go ; $j++){
	my $p = int(-1 * $temp[$j]->{pv}) ;
	print O "\t", $temp[$j]->{go}, " (p<1e-", $p, ")" ;
    }
    print O "\n" ;
}
