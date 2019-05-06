#!/bin/sh

 git clone https://github.com/lanl/nut.git -b feat-add-gtest nut
 cd nut
 mkdir build
 cd build
 cmake ..
 make VERBOSE=on -j 4 2>&1 | tee make.out
 ./test/nut_unittests
