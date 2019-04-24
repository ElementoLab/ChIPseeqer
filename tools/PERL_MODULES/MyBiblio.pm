package MyBiblio;
use strict;

use Bio::Biblio;
use XML::DOM;
use XML::XQL;
use XML::XQL::DOM;
use Database;
use File::Copy;

#   mysql db create here :
#     CREATE TABLE MYBIBLIO (ID int not null auto_increment primary key, AUTHORS
#
#
sub new {
    my ($self) = {};
    $self->{XML} = undef;
    
    $self->{DB}  = Database->new;
    $self->{DB}->connect("MYBIBLIO");
    my %dum = ();
    $self->{DBS} = \%dum;
    bless $self;
    return $self;
}


#
#  add an article + authors to a DB
#
sub db_add {
    my ($self, $p) = @_;

    # get all records with FIRSTAUTHOR ID 
    my $a_ref = $self->db_getArticlesByFirstAuthorAndYear($self->{FIRSTAUTHOR}, $self->{YEAR});

    # how many first
    my $nb    = scalar(@$a_ref);


    # create a FIRSTAUTHORID
    my @letters = ("", "a", "b", "c", "d", "e", "f", "g", "h");

    $self->{FIRSTAUTHORID}     = $self->{FIRSTAUTHOR} . $self->{YEAR} . $letters[$nb];

    
    #print "INSERT $self->{FIRSTAUTHORID}\n";
    
    

    # create SQL
    my $sql = $self->toSQL;
    
    
    # store !
    $self->{DB}->query($sql);


    # get the ID of the article
    my $aid = $self->{DB}->lastInsertId;
    
    # add the DB
    $self->db_addLink($aid);

    
    
    
    
}



sub db_getArticlesByDB {
    
    my ($self, $db) = @_;
    
    #$self->{DB}->setVerbose(1);
    my $sql = "SELECT * from ARTICLE2DB, ARTICLES WHERE DB = '$db' AND ARTICLES.ID = ARTICLEID ORDER BY FIRSTAUTHORID";
    
    my @a_tmp = $self->{DB}->queryAllRecords($sql);

    return \@a_tmp;

}



sub db_getArticleByPMID {
    
    my ($self, $a) = @_;
    
    #$self->{DB}->setVerbose(1);
    my $sql = "SELECT ID from ARTICLES WHERE PMID = '$a'";
    
    my $h_tmp = $self->{DB}->queryOneRecord($sql);

    return $h_tmp;

}


#
# check the existence of a link between an article and a DB
#
sub db_existsLink {
    
    my ($self, $aid, $db) = @_;
    
    my $sql = "SELECT * FROM ARTICLE2DB WHERE ARTICLEID = '$aid' AND DB = '$db'";
    
    my @a_tmp = $self->{DB}->queryAllRecords($sql);
    
    return scalar(@a_tmp);
}


sub db_exists {
    
    my ($self, $a) = @_;
    my $h_ref = $self->db_getArticleByPMID($a);
    return $h_ref->{ID};

}



sub db_getArticlesByFirstAuthorAndYear {
    my ($self, $a, $y) = @_;
    
    #$self->{DB}->setVerbose(1);
    my $sql = "SELECT FIRSTAUTHORID from ARTICLES WHERE FIRSTAUTHOR = '$a' and YEAR = '$y'";

    my @a_tmp = $self->{DB}->queryAllRecords($sql);

    return \@a_tmp;
}


#
#  add a DB entry for an article
#
sub db_addLink {
    my ($self, $id, $db) = @_;

    my $thisdb = ($db?$db:$self->{MYDB});

    
    die "Please define a DB ..\n" if (!$thisdb);
    
    my $sql = "INSERT INTO ARTICLE2DB (ARTICLEID, DB) VALUES ('$id', '$thisdb')";
    $self->{DB}->query($sql);
} 


#
#
#
sub setDB {
    my ($self, $p) = @_;

    $self->{MYDB} = $p;
}

#
#  loads DBs from a file
#
sub loadDBFile {
    my ($self, $f) = @_;
    
    open IN, $f or die "Cannot open DB file $f\n";
    
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l;
	$self->{DBS}->{$a[0]} = $a[1];
    }
	
    close IN;

    
    
}

#
#  recreate a whole BibTex database from a DB
#
sub generateBibTex {
    
    my ($self) = @_;
    
    my @a = localtime;
    

    if (-e $self->{DBS}->{$self->{MYDB}}) {
	copy $self->{DBS}->{$self->{MYDB}}, $self->{DBS}->{$self->{MYDB}} . ".$a[3]-$a[4]-$a[5]";
    }

    open OUT, ">" . $self->{DBS}->{$self->{MYDB}};

    # get all articles for this DB
    my $a_ref = $self->db_getArticlesByDB($self->{MYDB});
    
    foreach my $h_ref (@$a_ref) {
	$self->{PMID}          = $h_ref->{PMID};
	$self->{YEAR}          = $h_ref->{YEAR};
	$self->{VOLUME}        = $h_ref->{VOLUME};
	$self->{PAGES}         = $h_ref->{PAGES};
	$self->{ARTICLETITLE}  = $h_ref->{TITLE};
	$self->{NUMBER}        = $h_ref->{ISSUE};
	$self->{FIRSTAUTHORID} = $h_ref->{FIRSTAUTHORID};
	$self->{AUTHORS}       = $h_ref->{FASTAUTHORS};
	
	print OUT $self->toBibTex;
	print OUT "\n\n";

	
	
    }
    
    close OUT;
    
    # make a backup of the former one
    


    # convert to BibTex, save to file

    
    # add the manual ones at the end
    
}


#
#  get an article by PMID
#
sub getArticleByPMID {
    my ($self, $p) = @_;

    $self->{PMID}  = $p;

    $self->{XML} = new Bio::Biblio->get_by_id($p);

    #print $self->{XML};

    return undef if (!$self->{XML});

    my $parser = new XML::DOM::Parser;
    my $doc = $parser->parse($self->{XML}); 

    my @a_tmp             =  $doc->xql("MedlineCitation/Article/Journal/JournalIssue/PubDate/Year");
    $self->{YEAR}         = $a_tmp[0]->getChildNodes->item(0)->getData;
        

    $self->{JOURNAL}      = $doc->getElementsByTagName("MedlineTA")   ->item(0)->getChildNodes->item(0)->getData;
    $self->{ARTICLETITLE} = $doc->getElementsByTagName("ArticleTitle")->item(0)->getChildNodes->item(0)->getData;
    $self->{VOLUME}       = $doc->getElementsByTagName("Volume")      ->item(0)->getChildNodes->item(0)->getData;
    
    my $otmp = $doc->getElementsByTagName("Issue")->item(0);
    if ($otmp) {
	$self->{NUMBER}     = $otmp->getChildNodes->item(0)->getData;
    }
    
    
    $self->{PAGES}        = $doc->getElementsByTagName("MedlinePgn")  ->item(0)->getChildNodes->item(0)->getData;
    my $otmp = $doc->getElementsByTagName("AbstractText")->item(0);
    if ($otmp) {
	$self->{ABSTRACT}     = $otmp->getChildNodes->item(0)->getData;
    }

    # remove the dot
    $self->{ARTICLETITLE} =~ s/\.$//;
    
    my @a_authors = ();
    $self->{FIRSTAUTHOR} = undef;
    
    my $cnt = 0;
    foreach my $a ($doc->getElementsByTagName("Author")) {
	
	
	my $i = $a->getElementsByTagName("Initials")->item(0)->getChildNodes->item(0)->getData;
	my @a = split //, $i;
	foreach $_ (@a) { $_ .= "."; }
	$i    = join "", @a;
	my $l = $a->getElementsByTagName("LastName")->item(0)->getChildNodes->item(0)->getData;
	
	if ($cnt == 0) { $self->{FIRSTAUTHORID} = $l; }

	push @a_authors, "$i $l";

	$cnt ++;
    }

    $self->{FIRSTAUTHORID} .= $self->{YEAR};
    
    $self->{AUTHORS} = join(" and ", @a_authors);
        
}

#
#  
#
sub toSQL {
    my ($self) = @_;

    my $sql = "INSERT INTO ARTICLES  (
FIRSTAUTHORID,
TITLE,
PAGES,
ISSUE,
VOLUME,
YEAR,
FIRSTAUTHOR,
FASTAUTHORS,
ABSTRACT,
PMID)
VALUES (
'$self->{FIRSTAUTHORID}',
'$self->{ARTICLETITLE}',
'$self->{PAGES}',
'$self->{NUMBER}',
'$self->{VOLUME}',
'$self->{YEAR}',
'$self->{FIRSTAUTHOR}',
'$self->{AUTHORS}',
'$self->{ABSTRACT}',
'$self->{PMID}'
)
";
 

    return $sql;

    

}

sub toBibTex {
    my ($self) = @_;
    
   
    my $out = "\@Article{$self->{FIRSTAUTHORID},
  author =       {$self->{AUTHORS}},
  title =        {{$self->{ARTICLETITLE}}},
  journal =      {$self->{JOURNAL}},
  year =         {$self->{YEAR}},
  OPTkey =       {},
  volume =       {$self->{VOLUME}},
  number =       {$self->{NUMBER}},
  pages =        {$self->{PAGES}},
  OPTmonth =     {},
  OPTnote =      {},
  OPTannote =    {PMID:$self->{PMID}}
}
";
    
    return $out;
}


1;

