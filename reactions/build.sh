#!/bin/bash

echo 'Running configure'

##CINCLUDE="/opt/local/include"
##CLIBS="/opt/local/lib"
SUNDIR="/Users/bbose/Desktop/pnl/validation+code/build/src"
home_dir=`pwd`

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$home_dir/lib:$SUNDIR/lib64"

## Edit this command out if you just want to recompile src files
./configure --prefix=`pwd` LDFLAGS="-L$SUNDIR/lib64"  CPPFLAGS="-I$SUNDIR/include" CFLAGS='-w' CXXFLAGS='-g -std=c++11' LIBS="-lgsl -lstdc++ -lsundials_cvode -lsundials_nvecserial"
##./configure --prefix=`pwd` LDFLAGS="-L$CLIBS -L$SUNDIR/lib64 -Wl,-rpath,$SUNDIR/lib64"  CPPFLAGS="-I$CINCLUDE -I./include -I$SUNDIR/include" CFLAGS='-w' CXXFLAGS='-g -std=c++11' LIBS='-lgsl -lstdc++ -lsundials_cvode -lsundials_nvecserial'

echo "Running make"
make
make install

echo "Building dynamic library"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$home_dir/lib:$SUNDIR/lib64"
#export DYLD_LIBRARY_PATH="$PWD/lib:$CLIBS:$SUNDIR/lib64:$DYLD_LIBRARY_PATH"

g++ -I./include -L./lib  -L"$SUNDIR" -lgsl -lgslcblas  -lstdc++ -lcopter  -lsundials_cvode -lsundials_nvecserial -fPIC -shared  -Wall src/cosmo_react.cpp -o lib/libreact.so
##g++ -I./include -I$CINCLUDE -I$SUNDIR/include -L./lib -L$CLIBS -L$SUNDIR/lib  -lgsl -lgslcblas  -lstdc++ -lcopter  -lsundials_cvode -lsundials_nvecserial -fPIC -shared -Wall src/cosmo_react.cpp -o lib/libreact.so

echo "done."
