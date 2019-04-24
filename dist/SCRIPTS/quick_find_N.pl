use strict;

open IN, $ARGV[0];

my $l = <IN>;
my ($ch) = $l =~ /\>(.+?)$/;

my $pos = 0;
my $prev = undef;
my $t_end = undef;
my $t_start = undef;

while (my $l = <IN>) {

  chomp $l;
  $l =~ s/[\r\n\ ]*//g;
  $l = uc($l);
  my @a = split //, $l;

  foreach my $c (@a) {
    #print "/" . uc($c) . "/\n";

    if ($c eq 'N') {
      
     # print "N";

      if ($prev eq 'N') {

	$t_end++;  # increase

      } else { # end prev = N

        # prev not N, this just started an N region
        # print "prev is '$prev'\n";
	if (defined($prev) && defined($t_start)) {
	  print "$ch\t$t_start\t$t_end\n";
	}
	
	$t_start = $pos;
	$t_end   = $pos;
      } # end prev not N     

    } # end if $c == N

    $prev = $c;
    $pos++;
  }
  #print "\n";

}
close IN;

#if ($prev eq 'N') {
if ($t_end > $t_start) { 
print "$ch\t$t_start\t$t_end\n";
}
