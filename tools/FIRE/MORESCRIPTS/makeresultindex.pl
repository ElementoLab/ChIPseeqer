#!/usr/bin/perl -w

use lib "$ENV{FIREDIR}/SCRIPTS";

if (!$ENV{FIREDIR} || $ENV{FIREDIR} eq '') {
  die "The FIREDIR environment variable is not set. I will now set to $pwd for you. Please just retype the same command. See tutorial for more information.\n";
  $ENV{'FIREDIR'} = $pwd;
  exit;
}

# makeresultindex.pl : processes results for FIRE web interface

# declarations

use Sets;
use Table;
use File::stat;
use File::Copy;
use strict;


my $command = 'export PATH=$PATH:/usr/local/bin:/usr/local/netpbm/bin';
`$command`;

my $queryid = shift;
my $title   = shift;
if (!$title) {
  die "Please provide a title\n";
}
if (!$queryid) { exit; }

# functions

sub filelines
{
#usage &filelines(file)
my $file = $_[0];

open (FILE, "$file");
my @filelines = <FILE>;
close (FILE);
return scalar(@filelines);
}

sub sevenmerlist
{
  #usage @array = &sevenmerlist(file)
  my $filename = $_[0];
  my %motifs = ();
  
  open (FILE, "$filename");
  my @filelines = <FILE>;
  close (FILE);
  foreach my $line (@filelines) {
    my @currow = split("\t", $line);
    $motifs{$currow[8]} = $currow[0];
    if($currow[1] == 1) {
      $motifs{$currow[8]} =~ s/T/U/g;
    }
  }
  return %motifs;
}


sub sevenmerlistorder
{
  #usage @array = &sevenmerlist(file)
  my $filename = $_[0];
  my @motifs;
  
  open (FILE, "$filename");
  my @filelines = <FILE>;
  close (FILE);
  
  foreach my $line (@filelines) {
    my @currow = split("\t", $line);
    push(@motifs, $currow[8]);
  }
  
  return @motifs;
}


sub parsemimatrix {
  my ($f) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($f);
  my $a_ref = $ta->getArray();
  
  my @a = ();
  foreach my $r (@$a_ref) {
    my $c1  = 0;
    if (($r->[6] ne "") && ($r->[6] ne "nan") && ($r->[6] <= 100)) { 
      $c1 = 1;
    }

    my $c2 = 0;
    if (($r->[8] ne "") && ($r->[8] ne "nan") && ($r->[8] <= 100)) { 
      $c2 = 1;
    }
    if ($c1 == 1) { # || ($c2 == 1)) {
      my @b = ($r->[0], $r->[1]);
      push @a, \@b;
    }
  }

  print scalar(@a) . " interactions *****************\n";

  return \@a;
}

#
#  &combmotifmapdd(\$queryid, \$filepath{'baseurl'}, \'DNA_RNA/', $a_ref_int_DNARNA, $h_ref_sum_DNARNA);
#
sub combmotifmapdd {

  my ($queryidscalar, $baseurlscalar, $dirscalar, $a_ref_int, $h_ref_sum) = @_;
  my $queryid = Sets::filename($$queryidscalar);
  my $baseurl = $$baseurlscalar;
  my $dir = $$dirscalar;

  return if (@$a_ref_int == 0);
  
  print RESULTINDEX '<dd>
<dl>
<dt>Combined motif maps for co-localizing motifs</dt>
<dd>
<table>
<tr> 
  <td align="center">Motif 1</td>
  <td align="center">Motif 2</td>
  <td align="center">&nbsp;&nbsp;Location&nbsp;&nbsp;</td>
  <td align="center" colspan="3">Maps</td>
</tr>
';
  
  my $zcount = 0;
  
  foreach my $r (@$a_ref_int) {
    
    my $ff   = $h_ref_sum->{$r->[0]}->[8] . "_" . $h_ref_sum->{$r->[1]}->[8];     
    my $eps  = $$queryidscalar . "\_FIRE/$dir/$queryid" . '.mimatrix_OUT/' . $ff . '.eps';
    next if (! -e $eps);
    
    print RESULTINDEX '<tr';

    if ($zcount % 2 == 1) { 
      print RESULTINDEX ' class="z"'; 
    }

    print RESULTINDEX '>';
    
    my $m1 = $r->[0];
    my $m2 = $r->[1];

    if ($h_ref_sum->{$r->[0]}->[1] eq '1') {
      $m1 =~ s/T/U/g;
      $m2 =~ s/T/U/g;            
    }

    print RESULTINDEX '<td align="center">' . $m1 . '</td>';
    print RESULTINDEX '<td align="center">' . $m2 . '</td>';
    print RESULTINDEX "<td align=\"center\">" . ($h_ref_sum->{$r->[0]}->[1] eq '1' ?"3'UTR":"5'")  . "</td>";
    print RESULTINDEX '<td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.mimatrix_OUT/' . $ff . '.pdf';
    print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a></td><td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.mimatrix_OUT/' . $ff . '.png';
    print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a></td><td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.mimatrix_OUT/' . $ff . '.eps';
    print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a></td>' . '</tr>' . "\n";
    $zcount++;
  }
  
  print RESULTINDEX '</table>
</dd>
</dl>
</dd>
';

}

sub motifmapdd
{
  my ($queryidscalar, $baseurlscalar, $dirscalar, $motifhash, $motiforderarray, $h_ref_sum) = @_;
  my %motifs = %$motifhash;
  my @orderedmotifs = @$motiforderarray;
  my $queryid = Sets::filename($$queryidscalar);
  my $baseurl = $$baseurlscalar;
  my $dir = $$dirscalar;
  
  print RESULTINDEX '<dd>
<dl>
<dt>Motif maps</dt>
<dd>
<table>
<tr> 
  <td align="center">Motif</td>
  <td align="center">Location</td>
  <td align="center">&nbsp;&nbsp;&nbsp;Position&nbsp;&nbsp;&nbsp;<br/>bias</td>
  <td align="center">Orientation<br/>bias</td>
  <td align="center" colspan="3">Maps</td>
</tr>
';
  
  my $zcount = 0;
  
  foreach my $orderedmotif (@orderedmotifs) {
    
    my $eps = $$queryidscalar . "\_FIRE/$dir/$queryid" . '.summary_OUT/' . $orderedmotif . '.eps';
    next if (! -e $eps);
    
    print RESULTINDEX '<tr';

    if ($zcount % 2 == 1) { 
      print RESULTINDEX ' class="z"'; 
    }

    print RESULTINDEX '>';
    
    my $m = $motifs{$orderedmotif}; $m =~ s/U/T/g;

   

    print RESULTINDEX '<td>' . $motifs{$orderedmotif} . '</td>';
    print RESULTINDEX "<td align=\"center\">" . ($h_ref_sum->{$m}->[1] eq '1' ?"3'UTR":"5'")  . "</td>";
    print RESULTINDEX "<td align=\"center\">" . ($h_ref_sum->{$m}->[9] eq '1' ?'Y':'')  . "</td>";
    print RESULTINDEX "<td align=\"center\">" . ($h_ref_sum->{$m}->[10] eq '1' ?'Y':'') . "</td>";
    print RESULTINDEX '<td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.summary_OUT/' . $orderedmotif . '.pdf';
    print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a></td><td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.summary_OUT/' . $orderedmotif . '.png';
    print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a></td><td> <a href="';
    print RESULTINDEX $baseurl . $dir . $queryid . '.summary_OUT/' . $orderedmotif . '.eps';
    print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a></td>' . '</tr>' . "\n";
    $zcount++;
  }
  
  print RESULTINDEX '</table>
</dd>
</dl>
</dd>
';
}

sub basicfigs
{
my ($queryidscalar, $baseurlscalar, $dirscalar) = @_;

my $queryid = Sets::filename($$queryidscalar);
my $baseurl = $$baseurlscalar;
my $dir     = $$dirscalar;

print RESULTINDEX '<dd><acronym title="Finding Informative Regulatory Elements">FIRE</acronym> p-value matrix <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.summary.pdf';
print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.summary.png';
print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.summary.eps';
print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a>
</dd>

<dd><acronym title="Finding Informative Regulatory Elements">FIRE</acronym> density matrix <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.densities.pdf';
print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.densities.png';
print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.densities.eps';
print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a>
</dd>

<dd><acronym title="Finding Informative Regulatory Elements">FIRE</acronym> motif interaction matrix <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.fullmimatrix.pdf';
print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.fullmimatrix.png';
print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a> <a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.fullmimatrix.eps';
print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a>
</dd>

';
#print "$$queryidscalar, $$baseurlscalar, $$dirscalar\n";

if (-e "$$queryidscalar\_FIRE/$dir/clusters.eps") {

  print RESULTINDEX '
<dd><acronym title="Finding Informative Regulatory Elements">FIRE</acronym> cluster centroids and over-represented motifs <a href="';
  print RESULTINDEX "$baseurl$dir" . "clusters.pdf";
  print RESULTINDEX '"><acronym title="Portable Document Format">PDF</acronym></a> <a href="';
  print RESULTINDEX "$baseurl$dir" . "clusters.png";
  print RESULTINDEX '"><acronym title="Portable Network Graphic">PNG</acronym></a> <a href="';
  print RESULTINDEX "$baseurl$dir" . "clusters.eps";
  print RESULTINDEX '"><abbr title="Encapsulated PostScript">EPS</abbr></a>
</dd>
';
  
}

}

sub basicfiles
{
my ($queryidscalar, $baseurlscalar, $dirscalar) = @_;
my $queryid = Sets::filename($$queryidscalar);
my $baseurl = $$baseurlscalar;
my $dir = $$dirscalar;

print RESULTINDEX '<dd><a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.GOmotifs';
print RESULTINDEX '">Best <acronym title="Gene Ontology">GO</acronym> enrichments for motifs</a></dd>
<dd><a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.GOmotifs.full';
print RESULTINDEX '">All <acronym title="Gene Ontology">GO</acronym> enrichments for motifs</a></dd>
<dd><a href="';
print RESULTINDEX $baseurl . $dir . $queryid . '.motifreport';
print RESULTINDEX '">All motif occurrences of motifs predicted by <acronym title="Finding Informative Regulatory Elements">FIRE</acronym>, sorted by putative functionality</a></dd>
</dl>';

}

# variable declarations
my %filepath = ();
my %results = ();

# filename declarations

$filepath{'mainquerydir'} = "$queryid\_FIRE";
$filepath{'qmetafile'} = $filepath{'mainquerydir'} . '/meta.txt';
$filepath{'firedir'} = '.';
$filepath{'fireqfile'} = $queryid;
$filepath{'fireresult'} = "$queryid\_FIRE";


$filepath{'fire_results_DNARNA'} = $filepath{'fireresult'} . '/DNA_RNA';
$filepath{'fire_results_DNA'}    = $filepath{'fireresult'} . '/DNA';
$filepath{'fire_results_RNA'}    = $filepath{'fireresult'} . '/RNA';
my $ff = Sets::filename($queryid);
$filepath{'DNARNA_summary'}      = $filepath{'fire_results_DNARNA'} . '/' . $ff . '.summary';
$filepath{'DNA_summary'}         = $filepath{'fire_results_DNA'} . '/' . $ff . '.summary';
$filepath{'RNA_summary'}         = $filepath{'fire_results_RNA'} . '/' . $ff . '.summary';
$filepath{'baseurl'}             = './';
$filepath{'resultindex'}         = $filepath{'fireresult'} . '/index.htm';

$filepath{'DNARNA_mimatrix'}      = $filepath{'fire_results_DNARNA'} . '/' . $ff . '.mimatrix';
$filepath{'DNA_mimatrix'}         = $filepath{'fire_results_DNA'} . '/' . $ff . '.mimatrix';
$filepath{'RNA_mimatrix'}         = $filepath{'fire_results_RNA'} . '/' . $ff . '.mimatrix';

while ((my $key, my $value) = each (%filepath)){
#print "$key\t$value\n";
$value =~ s/([\n]|)+([\r]|)+//g;
$filepath{$key} = $value;
}

## make result index page

$results{'DNARNA_summary_count'} = &filelines($filepath{'DNARNA_summary'});
$results{'DNA_summary_count'}    = &filelines($filepath{'DNA_summary'});
$results{'RNA_summary_count'}    = &filelines($filepath{'RNA_summary'});

my %DNARNA  = &sevenmerlist($filepath{'DNARNA_summary'});
my %DNA     = &sevenmerlist($filepath{'DNA_summary'});
my %RNA     = &sevenmerlist($filepath{'RNA_summary'});

my @oDNARNA = &sevenmerlistorder($filepath{'DNARNA_summary'});
my @oDNA    = &sevenmerlistorder($filepath{'DNA_summary'});
my @oRNA    = &sevenmerlistorder($filepath{'RNA_summary'});

#
#  load summary files
#
my $ta = Table->new;

$ta->loadFile($filepath{'DNARNA_summary'});
my $h_ref_sum_DNARNA = $ta->getIndex(0);

$ta->loadFile($filepath{'DNA_summary'});
my $h_ref_sum_DNA    = $ta->getIndex(0);

$ta->loadFile($filepath{'RNA_summary'});
my $h_ref_sum_RNA    = $ta->getIndex(0);

#
#  load interaction files
#
my $a_ref_int_DNARNA = [];
$a_ref_int_DNARNA    = parsemimatrix($filepath{'DNARNA_mimatrix'})  if (-e $filepath{'DNARNA_mimatrix'});
my $a_ref_int_DNA    = [];
$a_ref_int_DNA       = parsemimatrix($filepath{'DNA_mimatrix'})     if (-e $filepath{'DNA_mimatrix'}   );
my $a_ref_int_RNA    = [];
$a_ref_int_DNA       = parsemimatrix($filepath{'RNA_mimatrix'})     if (-e $filepath{'RNA_mimatrix'}   );


system("cp $ENV{FIREDIR}/MORESCRIPTS/screen.css " . $filepath{'fireresult'});

open (RESULTINDEX, ">$filepath{'resultindex'}");
print RESULTINDEX '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
<title>FIRE: ';


print RESULTINDEX $title . '</title>
<link rel="stylesheet" type="text/css" href="screen.css" media="screen" /> 
</head>

<body>

<div id="page">

<table width="100%">
<tr><td>
<h1><acronym title="Finding Informative Regulatory Elements">FIRE</acronym></h1>
</td>
<td align="right">
<!-- <a href="http://tavazoielab.princeton.edu/FIRE/">Back to FIRE home page</a> -->
</td>
</tr>
</table>

<div id="main" class="normal">

<h2>';

print RESULTINDEX $title;

print RESULTINDEX '</h2>

<h3>Summary</h3>

<p><acronym title="Finding Informative Regulatory Elements">FIRE</acronym> discovered ';

if ($results{'DNARNA_summary_count'} > 0) {

print RESULTINDEX $results{'DNARNA_summary_count'};

print RESULTINDEX ' motifs from the query data, consisting of ';

print RESULTINDEX $results{'DNA_summary_count'};

print RESULTINDEX " from DNA (5') and ";

print RESULTINDEX $results{'RNA_summary_count'};

print RESULTINDEX " from RNA (3'UTRs).</p>
";

#
#  DNA/RNA 
#
if($results{'DNA_summary_count'} > 0 || $results{'RNA_summary_count'} > 0) {

  print RESULTINDEX '<h3>Combined <abbr title="Deoxyribonucleic acid">DNA</abbr> and <abbr title="Ribonucleic acid">RNA</abbr> results</h3>

<dl>
<dt>Figures</dt>';

  &basicfigs(\$queryid, \$filepath{'baseurl'}, \'DNA_RNA/');
  
  &motifmapdd(\$queryid, \$filepath{'baseurl'}, \'DNA_RNA/', \%DNARNA, \@oDNARNA, $h_ref_sum_DNARNA);

  &combmotifmapdd(\$queryid, \$filepath{'baseurl'}, \'DNA_RNA/', $a_ref_int_DNARNA, $h_ref_sum_DNARNA);

  print RESULTINDEX '
</dl>

<h3><abbr title="Deoxyribonucleic acid">Additional information for DNA motifs</abbr></h3>

<dl>';

&basicfiles(\$queryid, \$filepath{'baseurl'}, \'DNA/');

print RESULTINDEX '

<h3><abbr title="Ribonucleic acid">Additional information for RNA motifs</abbr></h3>

<dl>';

&basicfiles(\$queryid, \$filepath{'baseurl'}, \'RNA/');
  
} else {
  
if($results{'DNA_summary_count'} > 0) {

print RESULTINDEX '<h3><abbr title="Deoxyribonucleic acid">DNA results</abbr></h3>

<dl>
<dt>Figures</dt>';

&basicfigs(\$queryid, \$filepath{'baseurl'}, \'DNA/');

&motifmapdd(\$queryid, \$filepath{'baseurl'}, \'DNA/', \%DNA, \@oDNA, $h_ref_sum_DNA);

&combmotifmapdd(\$queryid, \$filepath{'baseurl'}, \'DNA/', $a_ref_int_DNA, $h_ref_sum_DNA);


print RESULTINDEX '<dt>Additional information</dt>';

&basicfiles(\$queryid, \$filepath{'baseurl'}, \'DNA/');

}

if($results{'RNA_summary_count'} > 0) {

print RESULTINDEX '<h3><abbr title="Ribonucleic acid">RNA results</abbr></h3>

<dl>
<dt>Figures</dt>';

&basicfigs(\$queryid, \$filepath{'baseurl'}, \'RNA/');

&motifmapdd(\$queryid, \$filepath{'baseurl'}, \'RNA/', \%RNA, \@oRNA);

&combmotifmapdd(\$queryid, \$filepath{'baseurl'}, \'RNA/', $a_ref_int_RNA, $h_ref_sum_RNA);

print RESULTINDEX '<dt>Files</dt>';

&basicfiles(\$queryid, \$filepath{'baseurl'}, \'RNA/');

}

}

} else {

print RESULTINDEX 'no motifs.</p>
';

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
