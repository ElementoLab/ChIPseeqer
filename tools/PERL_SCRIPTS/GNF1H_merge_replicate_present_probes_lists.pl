#
#  what does this program do ?
#

use strict;


my %H   = ();
my $cnt = 0;
open IN, $ARGV[0];

my %SAVE = ();

while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
 
    #$a[0] =  lc($a[0]);
    #$a[0] =~ s/_//g;
   
    $H{ $a[0] } ++;
    
    $SAVE{ $a[0] } = $l;

}
close IN;


open LOG, ">log_GNF1_merge_presents";

foreach my $h (keys(%H)) {

    if ($H{$h} > 1) {
	
	my %PRESENT = ();

	open IN, $ARGV[0];
	
	while (my $l = <IN>) {
	    chomp $l;
	    my @a = split /\t/, $l, -1;
	    my $t = shift @a; #$t = lc($t);
	    if ($t eq $h) {
		
		foreach my $p (@a) {
		    
		    $PRESENT{ $p } ++;
		    
		}
		
	    }
	    
	}
	close IN;
	
	
	print "$h";
	foreach my $p (keys(%PRESENT)) {
	    if ($PRESENT{ $p } == $H{$h}) {
		print "\t$p";
	    } else {
		print LOG "$p is present only $PRESENT{$p}/$H{$h} in $h\n";
	    }
	}
	print "\n";
	
    } else {
	print "$SAVE{$h}\n";
    }
}


close LOG;
