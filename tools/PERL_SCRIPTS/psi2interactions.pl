#
#  HYbrigenics to tab 
#

use XML::DOM;
use LWP::Simple;
use strict;

#print "$ARGV[0]\n";

my $parser = new XML::DOM::Parser;
my $doc    = $parser->parsefile($ARGV[0]);



#
#  load interactors
#

open OUT, ">$ARGV[0].seq";
my %H = ();
foreach my $proteinInteractor ($doc->getElementsByTagName("proteinInteractor")) {

  my $shortLabel    = $proteinInteractor->getElementsByTagName("shortLabel")->item(0)->getChildNodes->item(0)->getData;
  my $fullName      = $proteinInteractor->getElementsByTagName("fullName")  ->item(0)->getChildNodes->item(0)->getData;

  my $sequence      = $proteinInteractor->getElementsByTagName("sequence")  ->item(0)->getChildNodes->item(0)->getData;

  my $id            = undef; 
  
  #print "fn=$fullName\n";

  foreach my $secondaryRef ($proteinInteractor->getElementsByTagName("secondaryRef")) {
    
    my $db = $secondaryRef->getAttributeNode ("db")->getValue;
    
    #print "db=$db\n";

    if ($db =~ /GadFly V3/) {
      
      #print "$db\n";
      
      $id = $secondaryRef->getAttributeNode ("id")->getValue; 
    }
  }

  print OUT ">$fullName $shortLabel $id\n$sequence\n\n";

  $H{ $fullName } = $fullName;

}

close OUT;


my $interactions      = $doc->getElementsByTagName("interaction");

my $n                 = $interactions -> getLength;
#print $n;

for (my $i = 0; $i < $n; $i++) {
    

  #my @PMID = ();
  #my $bibref      = $interactions->item($i)->getElementsByTagName("bibref")->item(0);
  
  #next if !defined($bibref);

  # foreach my $priref ($bibref->getElementsByTagName("primaryRef")) {
  #	my $db = $priref->getAttributeNode ("db")->getValue; 
  #my $id = $priref->getAttributeNode ("id")->getValue; 
  
  #	if ($db eq "PubMed") {
  #push @PMID, $id;
  #}

	

  #}
    
  my @INT = ();
  foreach my $p  ($interactions->item($i)->getElementsByTagName("proteinInteractorRef")) {
    
    my $ref = $p->getAttributeNode ("ref")->getValue; 
    #$ref =~ s/ID_//;
    
    push @INT, $ref;
  }

  #
  # get the confidence value
  #
  my $confidence = $interactions->item($i)->getElementsByTagName("confidence")->item(0)->getAttributeNode("value")->getValue;

  my $confidence_txt = "Two Hybrid";
  if ($confidence ne "D") {
    $confidence_txt .= " (High Confidence)";
  }
  
  my $np = scalar(@INT);
    
  for (my $j=0; $j<$np-1; $j++) {
    for (my $k=$j+1; $k<$np; $k++) {
      print "$H{$INT[$j]}\t$H{$INT[$k]}\t$confidence_txt\t15710747";
      #print "\t"; print join("/", @PMID); 
      print "\n";
    }
  }
  
  
  
}


$doc->dispose;



