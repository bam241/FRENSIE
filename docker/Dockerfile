FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y git libssl-dev python-pip python-dev swig libboost-all-dev wget
RUN pip install numpy 

RUN apt remove --purge --auto-remove cmake
RUN version=3.17; \
    build=3;\ 
    mkdir ~/temp; \
    cd ~/temp; \
    wget https://cmake.org/files/v$version/cmake-$version.$build.tar.gz; \
    tar -xzvf cmake-$version.$build.tar.gz;\
    cd cmake-$version.$build/;\
    ./bootstrap; \
    make -j$(nproc);\ 
    make install