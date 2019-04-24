BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;


my $co = "$home/PROGRAMS/MIMOTIFS/REFSEQ/KELLIS_ALIGNMENTS/hg17_up2k_dwn2k_aligned.nmbased.fa";
my $a_ref_op = Sets::getFiles("*.optim7$ARGV[0].rep");


die "toto\n" if (! -e $co );

my $d = $ARGV[1];
die "no directory for the .nodups files ..\n" if (!defined($d));

foreach my $f (@$a_ref_op) {

  my $ff = $f;
  $ff =~ s/\.optim7$ARGV[0]\.rep//g;
  my $fh = $ff;
  $ff .= ".seeds7$ARGV[0]";

  die "toto\n" if (! -e $f );
  die "tata\n" if (! -e $ff);

  my $todo = "perl $home/PROGRAMS/MIMOTIFS/create_web_report.pl $f $ff $co $d/$fh > $fh.html";

  print "$todo\n";
  
  system($todo);


}


my $todo2 = "perl $home/PROGRAMS/MIMOTIFS/create_web_report_index.pl > index.html";
system($todo2);


