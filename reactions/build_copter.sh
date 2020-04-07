#!/bin/bash

./configure --prefix=`pwd` LDFLAGS="${LDFLAGS}" CFLAGS="${CFLAGS} -w" CXXFLAGS="${CXXFLAGS} -g -std=c++11" LIBS="-lgsl -lstdc++ -lsundials_cvode -lsundials_nvecserial"

echo "Running make"
make
make install
