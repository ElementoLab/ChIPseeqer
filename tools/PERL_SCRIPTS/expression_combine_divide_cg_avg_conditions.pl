BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;

if (scalar(@ARGV) == 0) {
  die "perl ~/PERL_MODULES/SCRIPTS/expression_combine_divide_cg_avg_conditions.pl --exps= --cons= --name=\n";
}

GetOptions ('cons=s'     => \$cons,
            'exps=s'     => \$exps,
	    'name=s'     => \$name);



my @a  = split /\,/, $exps;
my @b  = split /\,/, $cons;
 
my $n1 = scalar(@a);
my $n2 = scalar(@b);

my $txt1 = join(" ", @a);
my $txt2 = join(" ", @b);

my $outf1 = "$name" . "_MERGED";
my $todo = "perl ~/PERL_MODULES/SCRIPTS/expression_merge_microarray_conditions.pl $txt1 $txt2 > $outf1";
system($todo);

my $outf2 = $outf1 . "_RATIOS";
my $todo = "perl ~/PERL_MODULES/SCRIPTS/expression_get_average_fold_change.pl $outf1 $n1 $n2 > $outf2";
system($todo);

my $outf3 = $outf2 . "_CGED";
my $todo = "perl ~/PERL_MODULES/SCRIPTS/expression_drosophila_transform_affy_table.pl $outf2 > $outf3";
system($todo);

my $outf4 = $outf3 . "_AVGED";
my $todo = "perl ~/PERL_MODULES/SCRIPTS/expression_average_rows_with_same_id.pl $outf3 > $outf4.txt";
system($todo);


