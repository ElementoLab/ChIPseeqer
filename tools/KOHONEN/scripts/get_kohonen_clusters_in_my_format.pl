my $t = -1;
while (my $l = <STDIN>) {
 chomp $l;
 if ($l =~ /Cluster/) {
   $t ++ ;
 } else {
   print "$l\t$t\n"; 
 }

} 
