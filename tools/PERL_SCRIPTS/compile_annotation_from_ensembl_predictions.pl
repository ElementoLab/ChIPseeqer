use lib qw(/home/olly/PERL_MODULES);

use Sets;



my %h1 = ();
my %h2 = ();
my %h3 = ();

while (my $l = <STDIN>){
    chomp $l;
    
    my @a = split /\t/, $l;

    $a[6] =~ s/[\ \s\t]//g;
    $a[7] =~ s/[\ \s\t]//g;

    push @{ $h3{$a[4]} }, $a[5] if (($a[5] ne "") && !Sets::in_array($a[5], @{ $h1{$a[4]} }));
    #push @{ $h1{$a[4]} }, $a[6] if (($a[6] ne "") && !Sets::in_array($a[6], @{ $h1{$a[4]} }));
    #push @{ $h2{$a[4]} }, $a[7] if (($a[7] ne "") && !Sets::in_array($a[7], @{ $h2{$a[4]} }));

}


my @k1s = keys(%h1);
my @k2s = keys(%h2);
my @k3s = keys(%h3);



#my $a_ref = Sets::getUnionSet(\@k1s, \@k2s);

my $a_ref = \@k3s;

foreach my $k (@$a_ref) {
    print "$k\t";

    print join(",", @{ $h3{$k} }); 
    print "\n";

    #print join(",", @{ $h1{$k} }); print "\t";
    
    #print join(",", @{ $h2{$k} }); 
    
    #push @{ $h1{$a[4]} }, $a[6] if (($a[6] ne "") && !Sets::in_array($a[6], @{ $h1{$a[4]} }));

    
    #print "\n";

}
