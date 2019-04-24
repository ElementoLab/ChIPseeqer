package GO;
use strict;

# use a pile


sub new {
    
    my $self  = {};
    
    # single nucleotide counts
    my %h1               = ();
    $self->{NODES}       = \%h1;

    my %h2               = ();
    $self->{NODES_INDEX} = \%h2;

    bless($self);           # but see below
    return $self;
    
}

# GO:0008150: BP
# GO:0003674: MF
# GO:0005575: CC
sub load_OBO_Ontology {

  my ($self, $file) = @_;

  open IN, $file or die "cannot open $file\n";
  
  my @nodes = ();

  while (my $l = <IN>) {
    
    chomp $l;
    
    if ($l =~ /^\[Term\]/) {
      
      my $id   = undef;
      my $name = undef;
      my @is_a = ();
      
      while (1) {
      
	$l = <IN>; chomp $l;

	last if ($l eq "");
      
	# id
	if ($l =~ /^id\: (.+?)$/) {
	  $id = $1;
	}
      
	# is_a
	if ($l =~ /^is\_a\: (.+?) \!/) {
	  push @is_a, $1;
	}

	if ($l =~ /^relationship\: part\_of (.+?) \!/) {
	  push @is_a, $1;
	}
      
	# name 
	if ($l =~ /^name\: (.+?)$/) {
	  $name = $1;
	}
	
	
      }
      
      #print "$id\t$name\t"; print scalar(@is_a); print "\n";

      my %tmp = ("GO" => $id, "PREV" => \@is_a, "NAME" => $name);

      $self->{NODES}->{$id} = \%tmp;
    }

  }

  
  
  

  close IN;
  
  
}


sub printNames {
  my ($self) = @_;
  
  foreach my $n (keys( %{ $self->{NODES} } )) {
    print "$n\t$self->{NODES}->{$n}->{NAME}\n";
  }
  
}

sub printNodes {

  my ($self) = @_;

  foreach my $n (keys( %{ $self->{NODES} } )) {
    my $a_ref = $self->getAllParents($n);
    print "$n\t"; print join("\t", @$a_ref); print "\n";
  }

}


sub loadOntology {

    my ($self, $file) = @_; 
    
    my @pile = ();
    my $top  = 0;
    
    my  $cur_sps = 0;
    
    my @nodes = ();
    my $prev_go = undef;
    
    open IN, $file or die "GO.pm: Cannot open Ontology $file ..\n";
    
    while (my $l = <IN>) {
	
	#print $l;
	chomp $l;

	next if ($l =~ /^\!/);    
	
	# we read a line, count the number of spaces
	$l  =~ /^(\ *)/;  #print "found " . length($1) . " spaces ...\n";    
	
	my $new_sps = length($1);
	
	
	#
	# divide the line into units
	#
	
	my @as = split //, $l;
	my @units = ();
	my $it = '';
	foreach $l (@as) {
	    if ($l =~ /[\%\<\$]{1}/) {
		push @units, $it if (length(trim($it)) > 0);
		$it = '';
	    }
	    $it .= $l;
	}
	push @units, $it if (length(trim($it)) > 0);
	
	#print "Got " . scalar(@units) . " units\n";
	#print join("\n", @units);
	#print "\n";
	
	
	#
	# first unit
	#
	
	
	my $u = shift @units;
	
	$u =~ /[\%\<\$]\s*(.+?)\;/;
	my $name = $1;
	
	#print "name=$name\n";
	
	#
	# a unit may have multiple GO ids (slim)
	# 
	my @go = ();
	while ($u =~ /GO\:(\d{7})/g) {
	    push @go, int($1); 
	}
	
	#
	#
	#print ">> $l, curr is $cur_sps, new is $new_sps\n";
	if ($new_sps > $cur_sps) {
	    # add current node to the stack
	    push @pile, $prev_go;
	    #print "father of @go is " . $pile[$#pile] . "\n";
	    $cur_sps = $new_sps;
	    
	    
	    
	} elsif ($new_sps == $cur_sps) {	
	    $cur_sps = $new_sps;
	    #print "father of @go is " . $pile[$#pile] . "\n";
	    
	    
	} elsif ($new_sps < $cur_sps) {	
	    
	    for (my $i=0; $i<($cur_sps - $new_sps); $i++) {
		pop @pile;
	    }
	    #print "father of @go is " . $pile[$#pile] . "\n";
	    $cur_sps = $new_sps;
	}

	
	#
	# create as many new nodes as required
	# 
	foreach my $n (@go) {
	    
	    if (!$nodes[$n]) {
		#print "CREATE a new node at positon $n named $name\n";
		my @pre = ($pile[$#pile]); $n = int($n);
		my %tmp = ("GO" => $n, "PREV" => \@pre, "NAME" => $name);
		$nodes[$n] = \%tmp;
		
	    } else {
		if (!in_array($pile[$#pile], @{$nodes[$n]->{PREV}})) {
		    push @{$nodes[$n]->{PREV}}, int($pile[$#pile]);
		}
		
	    
	    }
	}
	
	
	if (scalar(@units) != 0) {

	    foreach my $u (@units) {

		
		$u =~ /[\%\<\$]\s*(.+?)\;/;
		my $name = $1;
		
		#print "name=$name\n";
		
		
		my $n = "";
		$u =~ /GO\:(\d{7})/;
		$n = int($1); 
		
		if (!$nodes[$n]) {
		    #print "CREATE a new node at positon $n named $name\n";
		    my @pre = ();
		    my %tmp = ("GO" => $n, "PREV" => \@pre, "NAME" => $name);
		    $nodes[$n] = \%tmp;
		    
		} 
		
		#
		# link current node to that node
		#
		
		foreach my $nn (@go) {
		    print "VOID\n" if ($nn eq "");
		    if (!in_array($n, @{$nodes[$nn]->{PREV}})) {
			push @{$nodes[$nn]->{PREV}}, int($n);
		    }
		}
		
		#
		# end linking
		#
		
	    
		
	    }
	}

	
	$prev_go = int($go[0]);
    }
    
    foreach my $n (@nodes) {
	
	$self->{NODES_INDEX}->{$n->{GO}} = 1;

	next if (!$n);
	
	#print "$n->{GO}\t";
	
	#print join("/", @{$n->{PREV}});
	
	
	#print "\t$n->{NAME}\n";
    
    }
    
    
    #print $nodes[$nodes[7126]->{PREV}->[0]]->{NAME} . "\n";;
    
    $self->{NODES} = \@nodes;
}


sub get_all_categories {
  my ($self) = @_;
  
  my @k = keys(%{ $self->{NODES_INDEX} });

  return \@k;
  
}


sub category_exists {
    my ($self, $c) = @_;

    return $self->{NODES_INDEX}->{$c};
}

sub print {
    
    my ($self) = @_;
    
    foreach my $n (@{$self->{NODES}}) {
	next if (!$n);
	
	print "$n->{GO}\t";
	print join("/", @{$n->{PREV}});
	print "\t$n->{NAME}\n";
    
    }
    
    
}



sub getName {
    
    my ($self, $i) = @_;
        
    return $self->{NODES}->{$i}->{NAME};
}


sub getFathers {
    
    my ($self, $i) = @_;
        
    return $self->{NODES}->{$i}->{PREV};
}



sub getAllParents {

    my ($self, $i) = @_;
    
    my @a = ();
    
    $self->getParents($i, \@a);
    
    my @b = ();

    # inverted inex
    my %ix = ();
    foreach my $aa (@a) {
	push @b, $aa if (!$ix{$aa});
	$ix{$aa} = 1;
    }

    #print join("\n", @b) . "\n";
    
    return \@b;
}


#
#  recursive function
#
sub getParents {

    
  my ($self, $i, $r) = @_;
  
  return if (!$self->getFathers($i));
  
  my $a_ref_p  = $self->getFathers($i);

  #print "parents of $i: " . join("*", @$a_ref_p); print "\n";

  foreach my $p (@$a_ref_p) {
    next if (!$p);
    push @$r, $p;
    $self->getParents($p, $r);
  }
  

}




sub getParentVector {

    my ($self, $a_ref) = @_;
    
    my @all = ();
    
    foreach my $p (@$a_ref) {
	
	my $a_ref_par = $self->getAllParents($p);
	
	push @all, @$a_ref_par;
    }
    
    push @all, @$a_ref;

   
    return removeDuplicates(\@all);
}



sub getParentVectorBef {

    my ($self, $a_ref) = @_;
    
    my @all = ();
    
    foreach my $p (@$a_ref) {
	
	my $a_ref_par = $self->getAllParents($p);
	
	push @all, @$a_ref_par;
    }
    
    push @all, @$a_ref;

  
    return \@all;
}


sub trim {
    my ($s) = @_;

    my $t = $s;
    $t =~ s/\ //g;

    return $t;
}


sub removeDuplicates {
    my ($a_ref) = @_;
    
    my @b = ();
    
    # inverted inex
    my %ix = ();
    foreach my $aa (@$a_ref) {
	push @b, $aa if (!$ix{$aa});
	$ix{$aa} = 1;
    }

    #print join("\n", @b) . "\n";
    
    return \@b;
    
}


sub in_array() {
    my $val = shift(@_);
 
    foreach my $elem(@_) {
        if($val == $elem) {
            return 1;
        }
    }
    return 0;
}



1;
