#!/bin/sh

 git clone https://github.com/lanl/nut.git -b feat-incorporate-3d-mesh nut
 cd nut
 mkdir build
 cd build
 cmake ..
 make VERBOSE=on -j 4 2>&1 | tee make.out
 pwd
 ls -la
 ./test/nut_unittests
