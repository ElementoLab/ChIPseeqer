BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;


my %H1 = ();
my %H2 = ();
my $l = <STDIN>;
print $l;
while (my $l = <STDIN>) {
  
   chomp $l;
   
   my @a = split /\t/, $l, -1;
   
   # remove ID
   my $c = shift @a;
   
   my @d = ();
   foreach my $r (@a) {
       if (($r ne "") && ($r ne 'NaN') && ($r ne 'nan')) {
	   push @d, $r;
       }
   }
   
    
   # average 
   my $avg = Sets::average(\@d);  #print "AVG=$avg\t";
   my @b   = ();

   foreach my $r (@a) {
     if (($r ne "") && ($r ne 'NaN') && ($r ne 'nan')) {
       if ($ARGV[0]) {
	 $r = log($r / $avg)/log(2.0);
       } else {
	 $r = $r / $avg;
       }
       push @b, sprintf("%4.3f", $r);
     } else {
       $r = "";
       push @b, $r;
     }
      
     
   }
   
   # divide, log
   print $c . "\t" . join("\t", @b) . "\n";
   
}

   
