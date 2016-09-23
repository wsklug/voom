autoheader && 
aclocal && 
automake -a --foreign &&
autoconf &&
./configure &&
cd src/  &&
make 

