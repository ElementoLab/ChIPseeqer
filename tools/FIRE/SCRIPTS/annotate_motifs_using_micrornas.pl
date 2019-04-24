use lib "$ENV{FIREDIR}/SCRIPTS";

#
# we want to know if our motifs cover 7nt out of the 8 final ones in miRNAs
#


use Table;
use Sets;
use Fasta;
use Getopt::Long;
use strict;

my $namefile     = undef;
my $micrornas    = undef;
my $summaryfile  = undef;
my $col          = 0;
my $zmax         = 3.0;  # min
my $doz          = 1;
if (@ARGV == 0) {
  die "perl annotate_motifs_using_micrornas.pl --micrornas=FILE --summaryfile=FILE --namefile=FILE\n";
}

GetOptions ('micrornas=s'       => \$micrornas,
	    'summaryfile=s'     => \$summaryfile,
	    'doz=s'             => \$doz,
	    'namefile=s'        => \$namefile);


#
#  read in micrornas
#

my @MIRNAS  = ();
my $fa      = Fasta->new;
$fa->setFile($micrornas);
while (my $a_seq = $fa->nextSeq()) {
  my ($n, $s) = @$a_seq; 
  $n =~ s/\ .+$//g;  
  $s =~ s/t/u/g if ($s =~ /u/);
  my $ss = uc($s);
  $ss =~ s/U/T/g;
  $ss = Sets::getComplement($ss);
  my @a = ($n, $s, $ss);
  push @MIRNAS, \@a;
}

my $ta = Table->new;

$ta->loadFile($summaryfile);
my $a_ref_sum = $ta->getArray();
my $a_ref_mot = $ta->getColumn(0);



my $h_ref_names = {};
if (defined($namefile)) {
  $ta->loadFile($namefile);
  $h_ref_names = $ta->getIndexKV(0,1);
} else {
  foreach my $r (@$a_ref_mot) {
    $h_ref_names->{ $r } = "-";
  }
}


my $cnt       = 1;
foreach my $r (@$a_ref_sum) {
    
  next if ($r->[1] != 1);  # only RNA motifs dealt with here (for now)

  my $re = $r->[$col];	   

  my $a_ref_out = &get_9mer_matches_mirna_library($re, \@MIRNAS);
  
  my $tn = scalar( @$a_ref_out );
  
  my @mirnas = ();
  foreach my $o (@$a_ref_out) {
    #print $o;
    my @a = split /\t/, $o, -1;
    $a[1] =~ s/^.+miR/miR/g;
    push @mirnas, $a[1];
  }

  if  ($tn > 0) {
  
   if ($doz == 1) {

    my @vals = ();
    for (my $i=0; $i<1000; $i++) {
      my $a_ref_shu = &get_shuffled_mirnas_library(\@MIRNAS);
      my $a_ref_sho = &get_9mer_matches_mirna_library($re, $a_ref_shu);
      
      my $n         = scalar( @$a_ref_sho );
      
      push @vals, $n;
    }
    
    my $avg = Sets::average(\@vals);
    my $std = Sets::stddev (\@vals);
    
    my $z   = ($tn - $avg) / $std;
    
    #print "z=$z\n";
    
    if ($z > $zmax) {
      $h_ref_names->{ $r->[0] } .= "/" . join("/", @mirnas);  
    }
   } else {
      $h_ref_names->{ $r->[0] } .= "/" . join("/", @mirnas);  
   }
  }
}
 

foreach my $r (@$a_ref_sum) {
  $h_ref_names->{$r->[0]} =~ s/^\-\///g;
  print "$r->[0]\t" . $h_ref_names->{$r->[0]} . "\n";
}

sub get_shuffled_mirnas_library {
  
  my ($a_ref_mirnas) = @_;
  
  my @aa = ();
  foreach my $m (@$a_ref_mirnas) {
    my $ss = Sets::shuffle_seq($m->[2]);
    my @a = ("", "", $ss);
    push @aa, \@a;
  }

  return \@aa;
}

sub get_9mer_matches_mirna_library {
  
  my ($re, $a_ref_mirnas) = @_;
 
  my $ire = $re;

  $re =~ s/^\.+//g;
  $re =~ s/\.+$//g;

  my $a_ref_re = Sets::get_array_from_re($re);

  my @CRE = ();
  push @CRE, $re;

  for (my $i=0; $i<@$a_ref_re-7; $i++) {
    my $myre = "";
    for (my $j=$i; $j<$i+7; $j++) {
      $myre .= $a_ref_re->[$j]; 
    }
    push @CRE, $myre;
  }

  my $cnt = 0;

  my @OUT = ();
  foreach my $myre (@CRE) {
    
    foreach my $m (@$a_ref_mirnas) {
      
      if ($cnt == 0) {
	
	if ($m->[2] =~ /($myre).{0,1}$/) {
	  
	  my $l    = length($1);
	  my $nn   = "N" x $l; 
	  my $mymi = $m->[2];
	  $mymi   =~ s/$myre/$nn/;
	  push @OUT, "$myre\t$m->[0]\t$ire\t$mymi\n";

	}

      } else {

	if ($m->[2] =~ /($myre)$/) {
	
	  my $l    = length($1);
	  my $nn   = "N" x $l; 
	  my $mymi = $m->[2];
	  $mymi   =~ s/$myre/$nn/;
	  push @OUT, "$myre\t$m->[0]\t$ire\t$mymi\n";
	  
	}
	
      }
      
    }
    $cnt ++;
  }

  return \@OUT;
  
}
