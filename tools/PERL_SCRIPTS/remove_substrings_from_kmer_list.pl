use lib qw(/home/olly/PERL_MODULES);
use Sets;

my $a_ref_kmer = Sets::readKmers($ARGV[0]);

my $cnt        = scalar(@$a_ref_kmer);

for (my $i=0; $i<$cnt-1; $i++) {
    
    
   for (my $j=$i+1; $j<$cnt; $j++) {

       next if (!defined($a_ref_kmer->[$j]));
       next if ((length($a_ref_kmer->[$j])-length($a_ref_kmer->[$j])) != 1); 
       my $pattern = $a_ref_kmer->[$j]->[0];
       my $cpattern = Sets::getComplement($pattern);
       
       #print "does " . $a_ref_kmer->[$i]->[0] . " contain " . $pattern . " ?\n";
       
       if (($a_ref_kmer->[$i]->[0] =~ /$pattern/) || ($a_ref_kmer->[$i]->[0] =~ /$cpattern/)) {
	   $a_ref_kmer->[$j] = undef;

	 print "does " . $a_ref_kmer->[$i]->[0] . " contain " . $pattern . " ? .. YES\n";
       }
    
   } 
    
}

foreach my $r (@$a_ref_kmer) {
    if ($r->[0]){
	print join("\t", @$r); print "\n";
    }
}
