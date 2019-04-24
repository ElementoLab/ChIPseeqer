#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

if (@ARGV == 0) {
  die "Args: gen ann pre dir\n";
}

my $gen = $ARGV[0];
my $ann = $ARGV[1];
my $pre = $ARGV[2];
my $dir = $ARGV[3];

my $masgen = Sets::insert_suffix_before_file_extension($gen, "_masked", ".fasta");

my $todo   = undef;

#$todo = "perl ~/PERL_SCRIPTS/fasta_mask_genes.pl $gen $ann > $masgen"; 
#print "$todo\n";
#system($todo) == 0 or die "Cannot exec: $todo\n";

$todo = "perl ~/PERL_SCRIPTS/extract_upstream_sequences_from_genome.pl --genome=$gen --annotation=$ann --noblast=1 --checkmaxlen=1 --lengthU=2000 --lengthD=0      > $dir/$pre\_u_2000_0.fa";
print "$todo\n";
system($todo) == 0 or die "Cannot exec: $todo\n";

$todo = "perl ~/PERL_SCRIPTS/extract_upstream_sequences_from_genome.pl --genome=$gen --annotation=$ann --noblast=1 --checkmaxlen=1 --lengthU=1000 --lengthD=0      > $dir/$pre\_u_1000_0.fa";
print "$todo\n";
system($todo) == 0 or die "Cannot exec: $todo\n";

$todo = "perl ~/PERL_SCRIPTS/extract_upstream_sequences_from_genome.pl --genome=$gen --annotation=$ann --noblast=1 --checkmaxlen=1 --lengthU=10000 --lengthD=10000 > $dir/$pre\_u_10000_10000.fa";
print "$todo\n";
system($todo) == 0 or die "Cannot exec: $todo\n";

$todo = "perl ~/PERL_SCRIPTS/extract_upstream_sequences_from_genome.pl --genome=$masgen --annotation=$ann --noblast=1 --checkmaxlen=1 --lengthU=1000 --lengthD=0   > $dir/$pre\_u_1000_0_masked.fa";
print "$todo\n";
system($todo) == 0 or die "Cannot exec: $todo\n";

$todo = "perl ~/PERL_SCRIPTS/extract_upstream_sequences_from_genome.pl --genome=$masgen --annotation=$ann --noblast=1 --checkmaxlen=1 --lengthU=2000 --lengthD=0   > $dir/$pre\_u_2000_0_masked.fa";
print "$todo\n";
system($todo) == 0 or die "Cannot exec: $todo\n";
