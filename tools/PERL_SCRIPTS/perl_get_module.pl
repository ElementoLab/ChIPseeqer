print <<EOF
package $ARGV[0];

sub new {
    my (\$self) = {};
    
    \$self->{XXXX}        = undef;
    \$self->{UNAME}       = `uname`; \$self->{UNAME} =~ s/\\n//g;
    bless(\$self);
    return \$self;
}

sub set {
    my (\$self, \$n)       = \@_;
    \$self->{XXXX}        = undef;
}

sub get {
    my (\$self)           = \@_;
    return \$self->{XXXX};
}





1;

EOF
