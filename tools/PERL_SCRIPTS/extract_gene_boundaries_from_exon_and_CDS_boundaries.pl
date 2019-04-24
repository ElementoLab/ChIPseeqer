use lib qw(/home/elemento/PERL_MODULES);

use Table;
use strict;
use Sets;

if (scalar(@ARGV) == 0) {
    die "usage : txt file must be : gene chr coding_exon_start coding_exon_end strand transcript_start transcript_end\n";
}

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();

my %H;
foreach my $r (@$a_ref) {
      
  if ($H{$r->[0]}) {
    
    #
    #  coding exon boundaries
    #
    if (($r->[2] ne "") && ($r->[3] ne "")) {
      $H{$r->[0]}->[2] = ( $H{$r->[0]}->[2] ne "" ? Sets::min($H{$r->[0]}->[2], $r->[2]) : $r->[2] );
      $H{$r->[0]}->[3] = ( $H{$r->[0]}->[3] ne "" ? Sets::max($H{$r->[0]}->[3], $r->[3]) : $r->[3] );
    }
    
    #
    #  transcript boundaries
    #
    if (($r->[5] ne "") && ($r->[6] ne "")) {
      $H{$r->[0]}->[5] = ( $H{$r->[0]}->[5] ne "" ? Sets::min($H{$r->[0]}->[5], $r->[5]) : $r->[5] );
      $H{$r->[0]}->[6] = ( $H{$r->[0]}->[6] ne "" ? Sets::max($H{$r->[0]}->[6], $r->[6]) : $r->[6] );
    }
    
  } else {

    my @a = ($r->[0], $r->[1], $r->[2], $r->[3], $r->[4], 
	     $r->[5], $r->[6]);
    $H{$r->[0]} = \@a;

    if (defined($r->[7])) {
	push @{ $H{$r->[0]} }, $r->[7];
    }

  }
  
}

foreach my $k (sort(keys(%H))) {
   
    # change exons transcripts into CDS
    my $r = $H{ $k };

 
    #my $t = $r->[5];
    #$r->[5] = $r->[2];
    #$r->[2] = $t;

    #my $t = $r->[6];
    #$r->[6] = $r->[3];
    #$r->[3] = $r->[6];

    if (($r->[2] eq "") || ($r->[3] eq ""))  {
      $r->[2] = $r->[5];
      $r->[3] = $r->[6];
    }

    print join("\t", @$r );
    print "\n";
}


