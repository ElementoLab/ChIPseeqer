#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

use Getopt::Long;

my $table = undef;
my $col   = 0;
my $NM    = undef;

GetOptions("table=s" => \$table,
           "col=s"   => \$col,
	   "NM=s"    => \$NM);
	   

my $ta = Table->new;
$ta->loadFile("$ENV{HOME}/PROGRAMS/PAGE/PAGE_DATA/ANNOTATIONS/human_go_orf_old/human_go_orf_old_genedesc.txt");
my $h_ref = undef;
if (!defined($NM)) { 
  $h_ref = $ta->getIndex(0);
} else {
  $h_ref = $ta->getIndex(1);
}

$ta->loadFile($table);
my $a_ref = $ta->getArray;


foreach my $r (@$a_ref) {
  my $cc = $r->[$col];  
  if (!defined($NM)) {
    $cc =~ s/^\ +//;
    print join("\t", @$r) . "\t" . $h_ref->{$cc}->[2] . "\n"; 
  } else {
    $cc =~ s/\-\d+$//;
    print join("\t", @$r) . "\t" . $h_ref->{$cc}->[0] . "\t" .  $h_ref->{$cc}->[2]  . "\n";
  }
}

