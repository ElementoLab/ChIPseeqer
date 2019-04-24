use strict;

#
#  read the first line
#

open IN, $ARGV[0];
my $l = <IN>; chomp $l;

my @b = split /\t/, $l, -1;
foreach my $r (@b) {
    next if ($r eq "");
    #print "$r\t";
    $r =~ s/\ //g;
    $r =~ s/^.+?\_//g;
    $r =~ s/\_[^\_]+$//g;
    
    #my @c = split /\_/, $r;
    #$r = $c[1];

    #print "$r\n";
}
    


my $n = scalar(@b);

my %PRESENT = ();

my $cnt = 0;
while (my $l = <IN>) {
     chomp $l;
     my @a = split /\t/, $l, -1;

     for (my $i=1; $i<$n; $i+=2) {
	 #print $b[$i] . "\t" . $b[$i+1] . "\n";

	 #print $a[$i] . "\t" . $a[$i+1] . "\n";
	 if ($a[$i+1] eq "P") {
	     
	     #die "oops i=$i, b[i+1]=." if ($b[$i+1] eq "");

	     push @{ $PRESENT{ $i+1 } }, $a[0];
	 }

     }

     $cnt ++;

     #last if ($cnt == 10);
     
}


close IN;


foreach my $k (sort(keys(%PRESENT))) {
    print $b[ $k ]; print "\t";  print join("\t", @{ $PRESENT{$k}}); print "\n";
}
