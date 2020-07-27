FROM ubuntu:18.04

RUN apt-get update; \
    apt-get install -y \
        cmake \
        vim \
        git \
        wget \
        doxygen \
        libpcre3 libpcre3-dev \
        gfortran \
        libssl-dev \
        python-pip \
        python-dev \
        libboost-all-dev \
        swig; \
    pip install numpy 

# get frensie
WORKDIR /root
COPY ./ /root/FRENSIE/
WORKDIR /root/FRENSIE
RUN mkdir bld
WORKDIR /root/FRENSIE/bld
RUN cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE ..
RUN make -j 8
