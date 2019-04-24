while ($l = <STDIN>) {
chomp $l;

@a = split /\ /, $l;



print "A\t" . join("\t", @a);
print "\n";
}
