FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y vim \
    git \
    autoconf \
    bison \
    flex \
     libpcre3 libpcre3-dev \
    gfortran \
    libssl-dev \
    python-pip \
    python-dev \
    wget
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
    make install;\
    cd .. ;\
    rm -rf cmake*




RUN git clone https://github.com/swig/swig.git;\
    cd swig/;\
    git fetch -all;\
    git checkout rel-4.0.0;\
    ./autogen.sh;\
    ./configure --prefix=/usr/local;\
    make -j8;\
    make install

RUN cd ~/temp;\
    wget https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz;\
    tar -xvf boost_1_72_0.tar.gz; \
    cd boost_1_72_0;\
    ./bootstrap.sh;\
    ./b2 install -s NO_BZIP2=1 link=shared runtime-link=shared cxxflags="-D_GLIBCXX_USE_CXX11_ABI=0" install;\
    cd ../ ;\
    rm -rf boost*

RUN apt-get install -y doxygen

# get frensie
WORKDIR /root
RUN echo "3"
RUN git clone https://github.com/bam241/FRENSIE.git
WORKDIR /root/FRENSIE
RUN git fetch
RUN git checkout docker
RUN mkdir bld
WORKDIR /root/FRENSIE/bld
run cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE ..
