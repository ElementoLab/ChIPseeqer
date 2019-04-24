#!/usr/bin/perl -w

# makeresultindex.pl : processes results for FIRE web interface

# declarations

use lib qw(/Genomics/sordune/www.8080/btsui/quantbio-tools/lib);
use strict;
use File::stat;
use File::Copy;

sub getFiles {
    my ($ls) = @_;



    my $s_list = `ls $ls`;

    #print "$s_list";

    my @a = split /\n/, $s_list;

    return \@a;


}


my $queryid = shift;
my $title   = shift;

open (RESULTINDEX, ">index.html");
print RESULTINDEX '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
<title>FIRE: ';


print RESULTINDEX $title . '</title>
<link rel="stylesheet" type="text/css" href="http://tavazoielab.princeton.edu/FIRE/screen.css" media="screen" /> 
</head>

<body>

<div id="page">

<table width="100%">
<tr><td>
<h1><acronym title="Finding Informative Regulatory Elements">FIRE</acronym></h1>
</td>
<td align="right"><a href="http://tavazoielab.princeton.edu/FIRE/">Back to FIRE home page</a></td>
</tr>
</table>

<div id="main" class="normal">

<h2>';

print RESULTINDEX $title;

print RESULTINDEX '</h2>

<h3>Summary</h3>
';

my $a_ref = &getFiles($queryid);
my @a = sort sortArrays @$a_ref;

foreach my $f (@a) {
  print RESULTINDEX '<p><a href="' . $f . '_FIRE/index.htm">' . $f . "</a></p>\n";
}


print RESULTINDEX '
<h3>Reference</h3>

<p>
<acronym title="Finding Informative Regulatory Elements">FIRE</acronym> is a motif discovery and characterization program based on mutual information.
</p>

</div>

<div id="footer">
<p>
<a href="http://tavazoielab.princeton.edu/">Tavazoie Lab</a>, <a href="http://www.genomics.princeton.edu/">Lewis-Sigler Institute for Integrative Genomics</a>, <a href="http://www.princeton.edu/">Princeton University</a>.
</p>
</div>

</div>
</body>

</html>';
close (RESULTINDEX); 


sub sortArrays {
  
  #my ($aa) = $a =~ /Arr\_(\d+)/;
  #my ($bb) = $b =~ /Arr\_(\d+)/;

  #my ($aa) = $a =~ /arr(\d+)/;
  #my ($bb) = $b =~ /arr(\d+)/;
  
  my ($aa) = $a =~ /stage(\d+)/;
  my ($bb) = $b =~ /stage(\d+)/;

  return $aa <=> $bb;
  #return $a cmp $b;
}
