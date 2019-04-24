#!/usr/bin/perl

use LWP::UserAgent;
use strict;


my $ua = new LWP::UserAgent;
use HTTP::Request::Common qw(POST);




my $req = $ua->request(POST("http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi",
			    [uid => $ARGV[0], dopt => 'MEDLINE',
			     cmd => 'Text', db => 'PubMed'
			     ])
		       );


my $content = $req->content();


#
# parse the content
#
$content =~ s/\<pre\>//g;
$content =~ s/\<\/pre\>//g;



my @TAGS = ();
my @CONT = ();
my $idx  = -1;
my @lines = split /\n/, $content;
foreach my $l (@lines) {
    
    # get the first chars

    my $ss1 = substr($l, 0, 6);
    my $ss2 = substr($l, 6);

    if ($ss1 =~ /^([A-Z]+)\ *\-\ /) {
	my $ctag = $1;
	$idx ++;
	$TAGS[$idx] = $ctag;

    }

    if ($CONT[$idx] ne "") {
	$CONT[$idx] .= " ";
    }
    $CONT[$idx] .= "$ss2";
}


my $title   = "";
my @authors = ();
my $volume  = "";
my $pages   = "";
my $journal = "";
my $year    = "";

my $n = scalar(@TAGS);

for (my $i=0; $i<$n; $i++) {
    if ($TAGS[$i] eq "TI") {
	$title = $CONT[$i];
	$title =~ s/\.$//;
    }

    if ($TAGS[$i] eq "PG") {
	$pages = $CONT[$i];
    }

    if ($TAGS[$i] eq "VI") {
	$volume = $CONT[$i];
    }

    if ($TAGS[$i] eq "TA") {
	$journal = $CONT[$i];
    }

    if ($TAGS[$i] eq "DP") {
	$year =  $CONT[$i];
	$year =~ s/\ .*$//g;
    }

    if ($TAGS[$i] eq "AU") {
	my $tmp = $CONT[$i];
	my ($t1, $t2) = $tmp =~ /^(.+)\ ([A-Z]+)$/;
	push @authors, "$t2 $t1";
    }
} 

my $firstauthorid = $authors[0];
$firstauthorid =~ s/^[A-Z]+\ //g;
$firstauthorid .= $year;

my $s_authors = join(" and ", @authors);

   
my $out = "\@Article{$firstauthorid,
  author =       {$s_authors},
  title =        {{$title}},
  journal =      {$journal},
  year =         {$year},
  OPTkey =       {},
  volume =       {$volume},
  OPTnumber =    {},
  pages =        {$pages},
  OPTmonth =     {},
  OPTnote =      {},
  OPTannote =    {PMID:$ARGV[0]}
}
";
    
print $out;

