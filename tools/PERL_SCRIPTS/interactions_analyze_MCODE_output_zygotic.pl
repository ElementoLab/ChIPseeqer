use lib qw(/home/elemento/PERL_MODULES);

open IN, $ARGV[0];

my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;


while (my $l = <IN>) {
  

  
  chomp $l;
  my @a = split /\t/, $l;
  print "CLUSTER $a[0]\n";
  
  my @b = split /, /, $a[4];

  open OUT, ">tmp";
  print OUT join("\n", @b);
  print OUT "\n";
  close OUT;

  #my $todo = "perl ~/PERL_MODULES/SCRIPTS/doGroupEnrichment.pl tmp ../../../DATA/DROSOPHILA/GO/go_full_index.txt ../../../DATA/DROSOPHILA/GO/GO_terms.txt -1";
  my $todo = "perl ~/PERL_MODULES/SCRIPTS/doGroupEnrichment.pl tmp ../mbt_index.txt -1 -1 0.001 1000000";
  system($todo);
  #<STDIN>;
}

close IN;
