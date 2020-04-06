#!/bin/bash

echo 'Running configure for COPTER'

./configure --prefix=`pwd` CFLAGS='-w' CXXFLAGS='-g -std=c++11' LIBS="-lgsl -lstdc++ -lsundials_cvode -lsundials_nvecserial"

echo "Running make"
make
make install
