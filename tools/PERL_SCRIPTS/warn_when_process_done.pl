while ( -e "/proc/$ARGV[0]") {
    sleep(60);
}
system("echo \"nothing\" | mail -s \"$ARGV[0] is done\" elemento\@princeton.edu");
