BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

my $file = undef;

while (my $l = <STDIN>) {
  if ($l =~ /^\!Sample_title = (.+?)$/) {
    
    my $file = $1;
    $file =~ s/[\,\'\"]//g;
    $file =~ s/[\/\ ]/\_/g;    
    open IN, ">$file";
  }

  if ($l !~ /^[\!\^\#]/) {
    print IN $l;
  }

}


close IN;
