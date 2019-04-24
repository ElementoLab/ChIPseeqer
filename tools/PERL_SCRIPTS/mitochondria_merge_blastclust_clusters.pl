BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use MyBlast;
use Fasta;

use strict;

if (@ARGV == 0) {
  die "usage : mitochondria_merge_blastclust_clusters.pl analysis+manual clusters all_seq.fa\n";
}

#
# read table of clusters to keep
#
my $a_ref_tokeep = Sets::readSet($ARGV[0]);
my @a_tokeep     = ();
my $cnt = 0;
foreach my $r (@$a_ref_tokeep) {
  if ($r =~ /^\*/) {
    $a_tokeep[$cnt] = 1;
  } else {
    $a_tokeep[$cnt] = 0;
  }

  $cnt ++;
} 


#
# read clusters, separate in 2 files 
#
my $ta = Table->new;
$ta->setDelim(" ");
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

my @a_goodclu = ();
my @a_badclu  = ();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($a_tokeep[$cnt] == 1) {
    push @a_goodclu, $r;
  } else {
    push @a_badclu , $r;
  }	
  $cnt ++;
}


#
#  create database of good clusters
#
my $create_database = 0;
my $dir             = "TMP";

my $se = Sequence->new;
$se->setBlastDB($ARGV[2]);    # all sequences

if ($create_database == 1) {
  if (! -e $dir) {
    mkdir $dir;
  }
  print "Outputing good clusters to $dir:\n";
  my $cnt = 0;
  foreach my $r (@a_goodclu) {
    
    my %TAKEN = ();
    open OUT, ">$dir/$cnt.fa";
    foreach my $g (@$r) {
      next if (defined($TAKEN{$g}));
      my $seq = $se->getSequenceFromBlastDB($g, 0, 0);
      if ($seq) {
	print OUT ">$g\n$seq\n\n";
	$TAKEN{ $g } = 1;
      } else {
	print "Cannot write $g.\n";
      }
    }
    close OUT;
    print "Cluster $cnt ... Done.\n";
    $cnt ++;
  }
}




my $mb = MyBlast->new;
$mb->setBlastProgram("blastp");
$mb->setFilter(0);
$mb->setEvalueThreshold("1");
$mb->setNbProcessors(2);
$mb->setExitOnCrash(1);


for (my $i=0; $i<@a_badclu; $i++) {

  my $r = $a_badclu[$i];
  
  #
  # traverse all genes in cluster
  #
  for (my $j=0; $j<@$r; $j++) {

    #
    # fetch gene, output to tmp file
    #
    my $g = $r->[$j]; next if ($g eq "");

    #next if ($g !~ /orf360\_Rhizopus\_oryzae/);
    #next if ($g !~ /atp8/);
    #next if ($g !~ /atp8_Metarhizium_an/);

    my $seq = $se->getSequenceFromBlastDB($g, 0, 0);
    next if (!$seq);

    my $tmpfile = Sets::getTempFile("toto");
    open OUT, ">$tmpfile" or die "cannot open $tmpfile\n";
    print OUT ">$g\n$seq\n";
    close OUT;
    
    $mb->setQueryDatabase($tmpfile);
    
    #
    # traverse all good clusters
    #
    
    my @THEHITS = ();

    for (my $k=0; $k<@a_goodclu; $k++) {

      #print "k=$k\n";

      # print "Blast $g to cluster $k\n";

      die "No $dir/$k.fa." if (! -e "$dir/$k.fa");
      
      $mb->setDatabaseDatabase("$dir/$k.fa"); 

      my $a_ref_hits = $mb->blastallMultiple;
      my $query_length = $mb->getQueryLength();

      if (@$a_ref_hits > 0) {
	#print "$g got " . scalar(@$a_ref_hits) . " hits to cluster $k.\n";
	my $good_hits = 0;
	my %species   = ();

	foreach my $hit (@$a_ref_hits) {
	  
	  my $name = $hit->{HIT_NAME};
	  my ($sp) = $name  =~ /\_([^\_]+\_[^\_]+)$/;
	  #print "name=$name, sp=$sp\n";
	  next if (defined($species{$sp}));

	  $species{$sp} = 1;
	  
	  my $a_ref_hsps = $hit->{HSPS};
	  my $hit_length = $hit->{HIT_LENGTH};
	  $a_ref_hsps = $mb->retain_non_overlapping_blocks($a_ref_hsps, 3);
	  my $a_ref_il = $mb->get_indentity_and_aligned_length($a_ref_hsps);
	  my ($I, $L) = @$a_ref_il;
	  #print "$I\t$L\n";
	  #print "maxl = " . Sets::max($query_length, $hit_length) . "\n";
	  my $LMIN = (2.0 / 3.0) * Sets::max($query_length, $hit_length);
	  
	  if (($L >= $LMIN) && ($I > 0.25))  {
	    $good_hits++;
	  }
	}
	#print "Got $good_hits good hits.\n";

	my @a_tmp = ($good_hits, $k);
	push @THEHITS, \@a_tmp;

      }

    }
    
    #
    # order the hits
    #
    @THEHITS = sort { $b->[0] <=> $a->[0] } @THEHITS;

    my $a_ref_bestclu = shift @THEHITS;
    
    if ($a_ref_bestclu->[0] > 1) {
      print "$g\t$a_ref_bestclu->[1]\t$a_ref_bestclu->[0]\n";
    } else {
      print "$g\t-\t-\n";
    }

	
    #
    #  
    #
    
    #
    #  order all HSPs, remove overlapping ones
    #
    #@$a_ref_hsps = sort { $a->{QFROM} <=> $b->{QFROM} } @$a_ref_hsps;
    #
	  
    #
    #  calculate the matching length, and the identity
    #
    #
    
    #
    
    #  # ge
	  
    #
    #
    #}
    #}
      
   
  
    unlink $tmpfile;
    #}
    
    #$cnt ++;
    ##last if ($cnt == 10);
    #}
    
    #if ($cnt_good > 0) {
    #  print "Compare to members of cluster $j ... ";
    #  print "$cnt_good matches / $cnt.\n";
    #}
    
  }
  
}

