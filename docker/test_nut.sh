#!/bin/sh

 echo "branch is "${branch}
 branch=$1
 git clone https://github.com/lanl/nut.git -b $branch nut
 cd nut
 git log -n 1 --oneline
 git branch
 mkdir build
 cd build
 cmake ..
 make VERBOSE=on -j 4 2>&1 | tee make.out
 pwd
 ls -la
 ./test/nut_unittests
