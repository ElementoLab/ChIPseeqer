cd PROGRAMS/pcre-7.4
./configure --prefix=$(pwd)/../../modules
make
make install
echo ""
echo "pcre-7.4 should now be installed in modules/lib/."
echo "To be sure, you can run './configure' again, and make sure you don't get any error messages,"
echo "or you can proceed to the compilation with 'make'."