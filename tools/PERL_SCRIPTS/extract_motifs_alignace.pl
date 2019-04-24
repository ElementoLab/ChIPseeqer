use lib qw(/home/olly/PERL_MODULES);



$f = $ARGV[0];
$n = $ARGV[1];


open IN, $f or die "toto\n";

@lines = <IN>;
$file  = join "", @lines;



my @a_parts = split /Motif \d*\n/, $file;



my @a_tmpl   = split /\n/, $a_parts[0];
foreach my $l (@a_tmpl) {
    if ($l =~ /^\#(\d+)\s+?(.+?)\s$/) {
	$index[$1] = $2;
    }
    
}


#print scalar(@a_parts) . "\n";;


#print "a=" . $a_parts[$n];

if (scalar(@a_parts) == 2) {
    @a_tmp1   = split /\n/, $a_parts[$n];
    #print $a_parts[$n] . "\n*****************************\n";
    
    
} else {
    
    @a_tmp1   = split /\n/, $a_parts[$n];
    #print $a_parts[$n] . "\n*****************************\n";
    
}

my $s_mot    = "Motif 1\n";


#print @a_tmp1;
foreach $l (@a_tmp1) {
    
    last if ($l =~ /MAP/);
    
    @a_tmp2 = split /\t/, $l;
    $s_first  = $a_tmp2[0];
    $s_secon  = $index[$a_tmp2[1]];
    #print "f=$s_first\n";


    $s_mot .= "$s_first" if ($s_first =~ /[ATCG]+/);

    foreach my $f (@$fs) {
	#next if $f->{DESCRIPTION} eq "transcription";
	
	 #$s_mot .= "$f->{DESCRIPTION}/";

    }

    
    #$s_mot  = "$s_first\n$s_mot" if  ($s_first =~ /\*/);
	
    if ($l =~ /\*/) {
	 $s_mot .= "$l";
	}
   
    $s_mot .= "\n";
 
    last if ($l =~ /MAP/);
    
}

$f =~ s/ace/motif/;


print "$s_mot";


