#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;
use Sets;
use strict;
use Table;

my $sets   = undef;
my $dbname = undef;
my $genelist = undef;

if (@ARGV == 0) {
  die "Args: --sets=STR --genelist=FILE --dbname=STR\n";
}

GetOptions("sets=s"     => \$sets,
           "genelist=s" => \$genelist,
	   "dbname=s"   => \$dbname);


my $ta = Table->new;
$ta->loadFile($genelist);
my $a_ref_genes = $ta->getColumn(0);



my $a_ref_files = Sets::getFiles($sets);

my @set_names = ();
my %INDEX = ();
foreach my $f (@$a_ref_files) {

  next if (($f eq "$dbname\_index.txt") || ($f eq "$dbname\_names.txt"));

  my $ff = $f;
  $ff =~ s/\.txt$//;

  push @set_names, $ff;
  
  my $h_ref_set = Sets::getIndex($f);

  foreach my $s (@$a_ref_genes) {
    push @{ $INDEX{ $s } }, $ff if (defined($h_ref_set->{$s}));
  }
  
}

my @gene_list      = @$a_ref_genes;
my @genes_in_index = keys(%INDEX);

my $a_ref_u   = Sets::getUnionSet(\@gene_list, \@genes_in_index);

die "Please define dbname\n" if (!defined($dbname));
 
open OUT1, ">$dbname\_index.txt";
open OUT2, ">$dbname\_names.txt";

foreach my $a (@$a_ref_u) {

  print OUT1 "$a\t";
  if (defined($INDEX{$a})) {
    print OUT1 join("\t", @{$INDEX{$a}});
  }
  print OUT1 "\n";

  
}
close OUT1;

foreach my $n (@set_names) {
  print OUT2 "$n\t$n\tP\n";
}
close OUT2;
