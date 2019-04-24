#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Getopt::Long;
use strict;


my $undef  = undef;
my $file   = undef;
my $dict   = undef;
my $k      = undef;
my $v      = undef;
my $col    = undef;
my $remdot = undef;
my $header = 1;

GetOptions ('table=s'  => \$file,
	    'dict=s'   => \$dict,
	    'col=s'    => \$col,
	    'k=s'      => \$k,
	    'v=s'      => \$v,
            'header=s' => \$header,
	    'remdot=s' => \$remdot,
	    'undef=s'  => \$undef);



if (!$file) {
    die "Usage : tr.. --table=s --dict=s --col=s --k=s --v=s --undef=s\n";
}


my $ta = Table->new;
$ta->loadFile($file);
my $a_ref = $ta->getArray();

$ta->loadFile($dict);
my $h_ref_kv = $ta->getIndexShifted();

#foreach my $g (keys(%$h_ref_kv)) {
#print "$g -> $h_ref_kv->{$g}\n";
#}

if ($header == 1) {
  my $a = shift @$a_ref;
  print join("\t", @$a); print "\n";
}


my %H = ();
foreach my $r (@$a_ref) {
  
  if (defined($h_ref_kv->{ $r->[ $col ] })) {
    
    my @ORFS = @{ $h_ref_kv->{ $r->[ $col ] } }; 
   
    my $n = shift @$r;

    foreach my $o (@ORFS) {
   
      $o =~ s/\ //g;
      next if ($o eq "");
      #if (!defined($H{$o}{$r->[1]})) { 
	print "$o\t" . join("\t", @$r) . "\n";
	#$H{$o}{$r->[1]} = 1;
      #}	

    }

  } else {
    #print STDERR "Lookup for $r->[$col] failed.\n";
  }

}



