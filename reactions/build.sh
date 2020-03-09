#!/bin/bash

echo 'Running configure'

sundial_dir="/Users/bbose/Desktop/pnl/validation+code/build/src"
home_dir=`pwd`

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$home_dir/lib:$sundial_dir/lib64"


## Edit this command out if you just want to recompile src files
#./configure --prefix=`pwd` LDFLAGS="-L$sundial_dir/lib64"  CPPFLAGS="-I$sundial_dir/include" CFLAGS='-w' CXXFLAGS='-g -std=c++11' LIBS="-lgsl -lstdc++ -lsundials_cvode -lsundials_nvecserial"

echo "Running make"
#make
#make install

echo "Building dynamic library"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$home_dir/lib:$sundial_dir/lib64"
g++ -I./include -L./lib  -L"$sundial_dir" -lgsl -lgslcblas  -lstdc++ -lcopter  -lsundials_cvode -lsundials_nvecserial -fPIC -shared  src/cosmo_react.cpp -o lib/libreact.so -Wall
