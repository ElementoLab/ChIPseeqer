#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;
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

#$ta->loadFile($file);
#my $a_ref = $ta->getArray();

open IN, $file or die "Cannot open $file.\n";

$ta->loadFile($dict);
my $h_ref_kv = $ta->getIndexKV($k, $v);

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


  my @b = split /\,/, $r->[$col];

  my @c = ();
  foreach my $bb (@b) {

    if ($h_ref_kv->{$bb} ne "") {
      push @c, $h_ref_kv->{$bb} if ( !Sets::in_array( $h_ref_kv->{$bb}, @c ) );
    }

  }

  if (@c > 0) {

    push @$r, join(",", @c);

    #if (defined($remdot)) {
    #  $r->[$col] =~ s/\.\d+$//;
    #}

    print join("\t", @$r); print "\n";

  } else {

        #print "$r->[$col]a .. pn\n";

        if ($undef eq "inf") {
	  push @$r, "inf";
	  print join("\t", @$r); print "\n";

        } elsif ($undef eq 'show') {
	  push @$r, $r->[$col];
	  print join("\t", @$r); print "\n";
	  
  	  #$r->[ $col ] = "inf"; 

        }
    }
    #print join("\t", @$r); print "\n";
}


close IN;
