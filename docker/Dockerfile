# Use an official Python runtime as a parent image
FROM tkzdockerid/coarct-testenv:clang-8-0-A

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
# ADD . /app

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get -y install googletest librandom123-dev

# RUN tar xvf clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz

# Define environment variable
ENV RANDOM123_DIR /usr
ENV BOOST_DIR /usr/include
ENV GTEST_DIR /usr/src/googletest/googletest

# End of file
