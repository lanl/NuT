## Run this down to === 1 === to make tkzdockerid/nut-testenv:clang-10-0-B-base
##FROM ubuntu:focal
##
### Set the working directory to /app
##WORKDIR /app
##
### Avoid chatty apt installs
##ARG DEBIAN_FRONTEND=noninteractive
##
### Copy the current directory contents into the container at /app
##ADD . /app
##
### Install any needed packages specified in requirements.txt
##RUN apt-get update && apt-get -y install cmake git clang-10 libtinfo5 \
##  googletest librandom123-dev \
##  docker.io vim xz-utils
##
### RUN cd /lib/x86_64-linux-gnu && ln -s libtinfo.so.5 libtinfo.so
##
### RUN tar xvf clang+llvm-10.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz
##
##RUN tar xvf gtest-master.tbz

## === 1 ===
## Run the rest to create tkzdockerid/nut-testenv:clang-10-0-B
FROM tkzdockerid/nut-testenv:clang-10-0-B-base

# Define environment variable
ENV RANDOM123_DIR /usr
ENV CC clang-10
ENV CXX clang++-10
ENV GTEST_SRC_DIR /app/master
ENV GTEST_ROOT /app/master/install
ENV PATH="/app/clang+llvm-10.0.0-x86_64-linux-gnu-ubuntu-18.04/bin:${PATH}"

# Now build googletest
# make and go to /app/master/build
RUN mkdir -p ${GTEST_SRC_DIR}/build && cd ${GTEST_SRC_DIR}/build && \
  cmake ${GTEST_SRC_DIR} -DCMAKE_INSTALL_PREFIX=$GTEST_ROOT -DCMAKE_BUILD_TYPE=Release && \
  make -j 4 install 2>&1 | tee make.out

# To build tkzdockerid/nut-testenv:clang-10-0-A
# FROM tkzdockerid/nut-testenv:clang-10-0-A-base


# End of file
