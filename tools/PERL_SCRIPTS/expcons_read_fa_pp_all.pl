BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Table;
use strict;

my $a_ref_chr = Sets::readSet($ARGV[0]);

my $ta = Table->new;

my $c = 1;

my $oo = Sets::get_array_from_re($ARGV[1]);
my $kk = @$oo;

my @out_seq = ();
my $out_pro = undef;
my $out_pos = undef;

my @flanks = (50, 125, 250, 500);

if (!defined($ARGV[2])) {

  $out_pro = "$ARGV[0].pro";
  $out_pos = "$ARGV[0].pos";
  foreach my $s (@flanks) {
    my $out_seq_f = "$ARGV[0]\_$s.seq";
    push @out_seq, $out_seq_f;
  }

} else {
    
  $out_pro = "$ARGV[2].pro";
  $out_pos = "$ARGV[2].pos";
  foreach my $s (@flanks) {
    my $out_seq_f = "$ARGV[2]\_$s.seq";
    push @out_seq, $out_seq_f;
  }

}

=pod

open OUTP, ">$out_pro" or die "cannot open file.\n";
open OUTO, ">$out_pos" or die "cannot open file.\n";

my @seq_h = ();
foreach my $f (@out_seq) {
  local *OUTS;
  open OUTS, ">$f" or die "cannot open file.\n";
  push @seq_h, *OUTS;
}

foreach my $chr (@$a_ref_chr) {
  
  # read conservation
  my $c_file = "$chr.pp";
  $ta->loadFile($c_file);
  my $h_ref = $ta->getIndexKV(0,1);

  # read alignment
  my $a_file = "$chr.fa.aln";
  my $fa = Fasta->new;
  $fa->setFile($a_file);
  
  my @ALN = ();
  my $seq = undef;
  my $cnt = 0;
  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;  
    my @a = split //, $s;
    
    if ($cnt == 0) {
      $seq = $s;
    }
    
    push @ALN, \@a;
    last;
    $cnt ++;
  }

  print "Read $chr.\n";

  $seq =~ s/\-//g;
  my @p_seq = split //, $seq;
  my $a_ref_pos = Sets::getREMotifPositions($ARGV[1], $seq);
  
  #my @a = @{ shift @ALN };

  my $a_a = $ALN[0];

  my $i = 0;
  my $j = 0;
  
  my $bt = undef;
  foreach my $r (@$a_a) {
    
    if ($r ne '-') {
      
      if (Sets::in_array($i, @$a_ref_pos)) {
	
	my $co = 0;
	my $st = "";
	
	for (my $k=$i; $k<$i+$kk; $k++) {
	  $co += $h_ref->{$k+1};
	  $st .= $p_seq[$k];
	}
	$co /= $kk;

	my $cnt_h = 0;
	foreach my $h (@seq_h) {
	  $bt = substr($seq, $i - $flanks[$cnt_h], 2 * $flanks[$cnt_h] + $kk);
	  print $h ">M$c\n$bt\n\n";
	  $cnt_h++;
	}	
	
	print OUTP "M$c\t" . sprintf("%3.2f", $co) . "\n";

	my $pp = $i + $kk;

	my $momo = substr($seq, $i, $i+$kk);
	my $como = Sets::getComplement($momo);
	my $strand = undef;
	if ($momo =~ /$ARGV[0]/) {
	  $strand =  1;
	} elsif ($como =~ /$ARGV[0]/) {
	  $strand = -1;
	}
	print OUTO "M$c\t$chr\t$i\t$pp\t$strand\n";
	$c ++;
      }
    
      
      
      $i++;
    }
    
    $j++;
  }    
  

}


foreach my $h (@seq_h) {
  close $h;
}
close OUTP;
close OUTO;

=cut

foreach my $f (@out_seq) {  
  my $todo = "formatdb -i $f -o T -p F";
  system($todo);
  $todo    = "perl ./detect_homologous_sequences.pl $f > $f.homologies";
  system($todo);  
}
