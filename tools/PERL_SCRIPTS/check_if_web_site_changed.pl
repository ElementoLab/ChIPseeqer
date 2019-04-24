use LWP::Simple;


my $h = get("http://www.amb-usa.fr/consul/niv/appointments/default.htm");

print $h;
