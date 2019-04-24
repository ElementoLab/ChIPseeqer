#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;
use Getopt::Long;
use strict;


my $undef  = undef;
my $file   = undef;
my $dict   = "$ENV{HOME}/PROGRAMS/ChIPseeqer/DATA/hg18/refGene.txt.07Jun2010.NM2ORF";
my $k      = 0;
my $v      = 1;
my $col    = 0;
my $remdot = undef;
my $header   = 1;
my $add      = 0;
my $multicol = 0;
my $idfile   = undef;

GetOptions ('table=s'    => \$file,
	    'dict=s'     => \$dict,
	    'col=s'      => \$col,
	    'k=s'        => \$k,
	    'v=s'        => \$v,
            'header=s'   => \$header,
	    'idfile=s'   => \$idfile,
	    'remdot=s'   => \$remdot,
	    'undef=s'    => \$undef,
	    'add=s'      => \$add,
	    'multicol=s' => \$multicol);



if (!$file) {
    die "Usage : tr.. --table=s --dict=s --col=s --k=s --v=s --undef=s\n";
}


my %OKID = ();
if (defined($idfile)) {
  open IN, $idfile or die "anoot open $idfile\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    $OKID{$a[0]} = 1;
  }
  close IN;


}

my $ta = Table->new;

#$ta->loadFile($file);
#my $a_ref = $ta->getArray();

open IN, $file or die "Cannot open $file.\n";

$ta->loadFile($dict);
my $h_ref_kv = {};
if ($multicol == 0) {
  $h_ref_kv = $ta->getIndexKV($k, $v);
} else {

  my $a_ref_di = $ta->getArray();
  foreach my $r (@$a_ref_di) {
    my $n = shift @$r;
          pop @$r while ($r->[ scalar(@$r) - 1 ] eq "");

    if (!defined($idfile)) {
      $h_ref_kv->{$n} = join(",", @$r);
    } else {
      my @newr = ();
      foreach my $s (@$r) {
	if (defined($OKID{$s})) {
	  push @newr, $s;
	}
      }
      if (@newr > 0) {
	$h_ref_kv->{$n} = join(",", @newr);
      } else {
	$h_ref_kv->{$n} = "N/A";
      }
    }
  }
}

#foreach my $g (keys(%$h_ref_kv)) {
#print "$g -> $h_ref_kv->{$g}\n";
#}

if ($header == 1) {
  #my $a = shift @$a_ref;
  #print join("\t", @$a); print "\n";
  my $l = <IN>;
  print $l;
}

while (my $l = <IN>) {
    
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $r = \@a;

  # split using ,
  my @b = split /\,/, $r->[$col];

  my @c = ();
  foreach my $bb (@b) {

    if ($h_ref_kv->{$bb} ne "") {
      push @c, $h_ref_kv->{$bb} if ( !Sets::in_array( $h_ref_kv->{$bb}, @c ) );
    }

  }

  if (@c > 0) {

    if ($add == 0) {
      $r->[ $col ] = join(",", @c);
    } else {
      $r->[ $col ] .= "\t" . join(",", @c);
    }
    #if (defined($remdot)) {
    #  $r->[$col] =~ s/\.\d+$//;
    #}

    print join("\t", @$r); print "\n";

  } else {

        #print "$r->[$col]a .. pn\n";

        if ($undef eq "inf") {
	  $r->[ $col ] = "inf";
	  print join("\t", @$r); print "\n";

        } elsif ($undef eq 'show') {
	  print join("\t", @$r); print "\n";
	  
  	  #$r->[ $col ] = "inf"; 

        }
    }
    #print join("\t", @$r); print "\n";
}


close IN;
