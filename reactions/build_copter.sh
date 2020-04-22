#!/bin/bash

./configure --prefix=`pwd` LDFLAGS="${LDFLAGS}" CFLAGS="${CFLAGS}" CXXFLAGS="${CXXFLAGS}"
echo "Running make"
make
make install
