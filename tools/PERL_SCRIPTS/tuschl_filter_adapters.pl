#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use bl2seq;
use strict;

my $bl = bl2seq->new;

my $fa = Fasta->new;
$fa->setFile($ARGV[2]);

my $ada1 = $ARGV[0];
my $fa1 = Fasta->new;
$fa1->setFile($ada1);
my $se1 = $fa1->nextSeq()->[1];

my $ada2 = $ARGV[1];
my $fa2 = Fasta->new;
$fa2->setFile($ada2);
my $se2 = $fa2->nextSeq()->[1];

my $tmpfile = Sets::getTempFile("/tmp/ada");

my $ii = 0;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;

  $s =~ s/\r//;
  open OUT, ">$tmpfile";
  print OUT ">$n\n$s\n";
  close OUT;

  #print "BEF: $s\n";
  #print "AD1: $se1\n";
  
  #print "AD2: $se2\n";
  if (0) {
  my $a_ref_m1 = $bl->bl2seq($ada1, $tmpfile);
  foreach my $r (@$a_ref_m1) {

    my $kill = 0;
    if ($r->[8] == 1) {
      #print "KILL1\n";
      $kill = 1;
    } else {
      #print "NO KILL1\n";
      $kill = 0;
    }

    if ($kill == 1) {
      my $s_st = 1; #$r->[8];
      my $s_en = $r->[9];
      my $ss   = substr($s, $s_st-1, $s_en - $s_st + 1);
      my $nn   = 'N' x length($ss);    
      substr($s, $s_st-1, $s_en - $s_st + 1)  = $nn;    
    }
    last;
  }
}

  if (1) {
    my $a_ref_m2 = $bl->bl2seq($ada2, $tmpfile);
    foreach my $r (@$a_ref_m2) {
      my $s_st = $r->[8];
      my $s_en = $r->[9];

      my $aa   = substr($s, $s_en);

      my $kill = 0;
      # kill if mostly As to the right
      if (nonanum($aa) <= 1) {
	
	#print "KILL2\n";
	$kill = 1;
      }
      # otherwise kill if big fragment (>=5) matching the begnning of adapter (<=2)
      elsif (($r->[6] <= 2) && ($r->[3] >= 5)) {
	#print "KILL2\n";
	$kill = 1;
      } 
      # if e<1e-4
      elsif ($r->[10] <= 1e-4) {
	#print "KILL2\n";
	$kill = 1;

      } else {
	#print "NO KILL2\n";
	$kill = 0;
      }	

      if ($kill == 1) {
	my $ss   = substr($s, $s_st-1);
	my $nn   = 'N' x length($ss);    
	substr($s, $s_st-1)  = $nn;    
      }
      last;
    }    
  }
  
  my $mys = $s;
  $mys =~ s/^N+//;
  $mys =~ s/N+$//;
  my $af = &afrac($mys);
  #print "$af\n";
  if ((length($mys) >= 20) && ($af < 1)) {
    print ">$n\n$mys\n";  
    print "\n";
  }

  $ii ++;

  if ($ii % 1000 == 0) {
    print STDERR "******************************** $ii reads.\n";
  }
  #<STDIN>;
}


sub afrac {
  my ($s) = @_;
  my @a = split //, $s;
  if (@a == 0) {
    return 0;
  }
  my $nn = 0;
  my $aa = 0;
  foreach my $nt (@a) {
    if ($nt eq 'A') {
      $aa ++;
    }
    if ($nt =~ /[ACGT]/) {
      $nn ++;
    }
  }
  return $aa/$nn;
}


sub nonanum {
  my ($s) = @_;
  my @a = split //, $s;
  if (@a == 0) {
    return 0;
  }
  my $aa = 0;
  foreach my $nt (@a) {
    if ($nt ne 'A') {
      $aa ++;
    }
  }
  return $aa;
}
