autoheader && 
aclocal && 
automake -a --foreign &&
autoconf &&
./configure &&
make 

