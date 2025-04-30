FROM ubuntu:24.04

ARG is_dev

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gcc-14 g++-14 \
    make \
    build-essential \
    cmake \
    git \
    zlib1g-dev
RUN update-alternatives --install /usr/bin/g++ g++ /bin/g++-14 14 
RUN update-alternatives --install /usr/bin/gcc gcc /bin/gcc-14 14

WORKDIR /
RUN if [[ -z "$is_dev" ]] ; then git clone -b dev https://github.com/rmcolq/charon.git; else git clone https://github.com/rmcolq/charon.git; fi
WORKDIR /charon/build

RUN cmake -DCMAKE_BUILD_TYPE=RELEASE .. > cmake.log 
RUN make -j4 > make.log
RUN make install > make_install.log

WORKDIR /

SHELL ["/bin/bash", "-c"]
