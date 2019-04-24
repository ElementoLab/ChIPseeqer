BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use XML::DOM;
use LWP::Simple;
use Log;
use Sets;
use strict;

#print "$ARGV[0]\n";

my $parser = new XML::DOM::Parser;
my $doc    = $parser->parsefile($ARGV[0]);


my $log = Log->new;
#
#  load interactions
#


#
#  traverse entries
#
foreach my $entry ($doc->getElementsByTagName("entry")) {
  
  #
  # get the PMID refs and exp desc
  #

  my @PMIDS = ();
  my @EXPES = (); 
  foreach my $experimentDescription ($entry->getElementsByTagName("experimentDescription")) {
    
    my $primaryRef = $experimentDescription->getElementsByTagName("primaryRef")->item(0);
  
    my $db = $primaryRef->getAttributeNode ("db")->getValue; 
    my $id = $primaryRef->getAttributeNode ("id")->getValue; 

    if ($db eq "PubMed") {
      push @PMIDS, $id if (!Sets::in_array($id, @PMIDS));
    }
    

    #
    # get the experiment description
    #

    my $method = $experimentDescription->getElementsByTagName("interactionDetectionMethod")->item(0)
      ->getElementsByTagName("names")->item(0)
	->getElementsByTagName("shortLabel")->item(0)
	  -> getChildNodes->item(0)->getData;
    
    #my $shortLabel = $experimentDescription->getElementsByTagName("shortLabel")->item(0)->getChildNodes->item(0)->getData;

    push @EXPES, $method if (defined($method));

  }

 

  

  my @TO = ();
  foreach my $proteinInteractor ($entry->getElementsByTagName("proteinInteractor")) {
    
    # get the label ..
    my $shortLabel = $proteinInteractor->getElementsByTagName("shortLabel")->item(0)->getChildNodes->item(0)->getData;
    
    # get the aliases ..
    my @aliases = ();
    foreach my $alias ($proteinInteractor->getElementsByTagName("alias")) {
      push @aliases, $alias->getChildNodes->item(0)->getData;
    }

    # get the Entrez Protein ID
    my $primaryRef = $proteinInteractor->getElementsByTagName("primaryRef")->item(0);
  
    my $db = $primaryRef->getAttributeNode ("db")->getValue; 
    my $id = $primaryRef->getAttributeNode ("id")->getValue; 

    if ($id) {
      my @a_tmp = ($shortLabel, $id, \@aliases);
      push @TO, \@a_tmp;
    }
    
  }
 
  if (scalar(@TO) == 2) {
    print "$TO[0]->[1]";

    #print "\t" . join("/", @{ $TO[0]->[2] }) . "\t"; 

    print "\t$TO[1]->[1]";

    #print "\t" . join("/", @{ $TO[1]->[2] }); 

    
    #if (scalar(@EXPES) > 0) {
    print "\t" . join("/", @EXPES);
    #}
    
    #print "\tTO CHECK";
    
    #if (scalar(@PMIDS) > 0) {
    print "\t" . join("/", @PMIDS);
    #}


    print "\n";

  } else {
    
    my $n = scalar(@TO);
    $log->log("$n interactions ..\n");
    
  }
  
}

$doc->dispose;



