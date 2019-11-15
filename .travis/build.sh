#!/bin/bash

mkdir build
cd build
mkdir w2rap-contigger
cmake -DCMAKE_INSTALL_PREFIX=w2rap-contigger .. ${CMAKE_OPTIONS}
make all
make install

tar cz w2rap-contigger > w2rap-contigger-${TRAVIS_OS_NAME}.tar.gz