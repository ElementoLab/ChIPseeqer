use strict ;

my ($expfile, $indexfile, $category, $cluster, $exptype) = @ARGV ;

my @g_in_c ;
my $min = $cluster ;
my $max = $cluster ;

if ($exptype eq "continuous"){
    $min = substr($cluster, index($cluster, "[")+1, index($cluster, " ")- index($cluster, "[") -1) ;
    $max = substr($cluster, index($cluster, " ")+1, index($cluster, "]")- index($cluster, " ") -1) ;
}

open I, "< $expfile" or die $expfile;
while (my $l = <I>){
    chomp $l;
    my ($g, $v) = split(/\t/, $l) ;
    if ($v <= $max and $v>=$min){
	push(@g_in_c, $g) ;
    }
}

my @genes ;
my (@des) = split(/\s/, $category) ;
my $cat = shift (@des) ;
if (! defined $cat){
    $cat = $category ;
}

open I, "< $indexfile" or die;
while (my $l = <I>){
    chomp $l;
    my ($g, @f) = split(/\t/, $l) ;
    next if (! (grep {$_ eq $cat} (@f))) ;
    next if (! (grep {$_ eq $g} (@g_in_c))) ;
    print $g, "\n" ;
}
