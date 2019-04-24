use lib qw(/home/olly/PERL_MODULES);
use Sets;

my $a_ref_kmer = Sets::readKmers($ARGV[0]);

my $cnt        = scalar(@$a_ref_kmer);

for (my $i=0; $i<$cnt-1; $i++) {
    
    
   for (my $j=$i+1; $j<$cnt; $j++) {
       
       next if (length($a_ref_kmer->[$i]->[0]) !=
		length($a_ref_kmer->[$j]->[0]));

       my @a = split //, $a_ref_kmer->[$i]->[0];
       my @b = split //, $a_ref_kmer->[$j]->[0];
       my @d = split //, Sets::getComplement($a_ref_kmer->[$j]->[0]);
       my @c = ();
       my $len = scalar(@a);
       
       my $score1 = undef;
       my $score2 = undef;
       

       #  first sense
       my $cnt = 0;
       for (my $k=0; $k<$len; $k++) {
	   if ($a[$k] ne $b[$k]) {
	       $c[$k] = "[" . $a[$k] . $b[$k] . "]";
	       $cnt ++;
	   } else {
	       $c[$k] = $a[$k];
	   }
       }
       
       my $re = join("", @c);

       if ($cnt <= 1) {
	   print "$a_ref_kmer->[$i]->[0] $a_ref_kmer->[$i]->[4]/ $a_ref_kmer->[$j]->[0] $a_ref_kmer->[$j]->[4] => $re ";
	   
	   
	   my $txt = `../FASTCOMPARE/recompare -re \"$re\" -fasta1 ../FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/CER.seq -fasta2 ../FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/BAY.seq -nbgenes 4358`;
	   ($score1) = $txt =~ /\=\ ([\d\.]+)/;
	   
	   print "$score1";

	   print "**"  if ($score1 > $a_ref_kmer->[$i]->[4]);
	   
	   print "\n";
       }

       #  second sense
       my $cnt = 0;
       for (my $k=0; $k<$len; $k++) {
	   if ($a[$k] ne $d[$k]) {
	       $c[$k] = "[" . $a[$k] . $d[$k] . "]";
	       $cnt ++;
	   } else {
	       $c[$k] = $a[$k];
	   }
       }
       
       my $re = join("", @c);

       if ($cnt <= 1) {
	   print "$a_ref_kmer->[$i]->[0] $a_ref_kmer->[$i]->[4]/ $a_ref_kmer->[$j]->[0] (C) $a_ref_kmer->[$j]->[4] => $re ";
	       
      
	   my $txt = `../FASTCOMPARE/recompare -re \"$re\" -fasta1 ../FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/CER.seq -fasta2 ../FASTCOMPARE/WEBSITE/RESULTS/CER_BAY/DATA/BAY.seq -nbgenes 4358`;
	   ($score2) = $txt =~ /\=\ ([\d\.]+)/;

	   print "$score2";

	   print "**"  if ($score2 > $a_ref_kmer->[$i]->[4]);

	   print "\n";
       }

       
       # remove the kmer that improved 
       if (($score1 > $a_ref_kmer->[$i]->[4]) || ($score2 > $a_ref_kmer->[$i]->[4])) {
	   
	   $a_ref_kmer->[$i] = undef;
	   
       }
       
       
          
   }
    

}
