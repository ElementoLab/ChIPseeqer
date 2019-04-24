BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use Table;
use strict;

my $mapfile       = Sets::get_parameter(\@ARGV, "-mapfile");
my $annotfile     = Sets::get_parameter(\@ARGV, "-annotfile");
my $minnbprobes   = Sets::get_parameter(\@ARGV, "-minnbprobes");
if (!$minnbprobes) {
  $minnbprobes = 8;
}

my $h_ref_TG = undef;
if ($annotfile) {
  my $ta = Table->new;
  $ta->loadFile($annotfile);
  $h_ref_TG = $ta->getIndexKV(1,0);
}

#
#  combine all cDNAs of a same probe set
#
my %H = ();
open IN, $mapfile;
while (my $l = <IN>) {
    chomp $l;    
    my @a = split /\t/, $l;
    $H{ $a[0] }->{ $a[1] } ++ if ($a[2] == 25);
}
close IN;



my @ps = keys(%H);
my %PROBE_SETS = ();
open LOG, ">log_affy_mapping.txt";

#
# traverse the probe set names
#
foreach my $p (@ps) {

  #
  #  get sets of mRNAs for this Transcript
  #
  my @ts = keys(%{ $H{$p} });
  my @ma = ();

  #
  #  for a given probe set, gather all CGs, with the number of probes matching the CG
  #
  foreach my $t (@ts) {
    my $n =  $H{ $p }->{ $t };
    my @tm = ($t, $n); 
    push @ma, \@tm;
  }
    
    #
    #  sort based on the number of probes
    #
    @ma = sort { $b->[1] <=> $a->[1] } @ma;
    my $nbt = scalar(@ma);


    #
    #  assign valid probe sets to genes, only one probe set per gene 
    #
    my %TAKEN = ();    
    my @txt = ();

    for (my $i=0; $i<$nbt; $i++) {
	
	my $t = $ma[$i][0]; 
	
	if (defined($annotfile)) {
	  if (defined($h_ref_TG->{ $t })) {
	    $t = $h_ref_TG->{ $t };
	  } else {
	    next;
	  }
	}

	my $n = $ma[$i][1];
	
	if (!defined($TAKEN{$t}) && ($n >= $minnbprobes) && !Sets::in_array($t, @txt)) {
	    push @txt, $t;
	    $TAKEN{ $t } = 1;
	}
    }
    
    if (scalar(@txt) > 0) {
	print "$p\t"; print join("/", @txt); print "\n";
    } else {
	print LOG "$p could not be assigned to any gene ..\n";
    }

}
close LOG;



