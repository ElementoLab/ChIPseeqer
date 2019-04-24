use lib qw(/home/elemento/PERL_MODULES);
use Sets;

open IN, $ARGV[0];
while (my $l = <IN>)  {
  
  if ($l =~ /^\*a (.+)$/) {

    if (scalar(@stack) > 0) {
      my $n = $stack[0];
      


      my $a_ref = Sets::removeDuplicates(\@stack);

      print "$n\t"; print join("\t", @$a_ref); print "\n";
    }
    
    @stack = ();

    #if ($l =~ /^\*a (CG\d+)/) {
    my $a = $1;
    if ($a !~ /D...\\/) {
      push @stack, $a;
    }
    #}
    
  }
  
  if ($l =~ /^\*i (CG\d+)/) {
    push @stack, $1;
  }

}
close IN;
