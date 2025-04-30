FROM ubuntu:24.04

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

RUN mkdir -p /charon
COPY src /charon/src
COPY include /charon/include
COPY cmake /charon/cmake
COPY lib /charon/lib
COPY .git /charon/.git
COPY CMakeLists.txt version.h.in /charon
WORKDIR /charon/build

RUN cmake -DCMAKE_BUILD_TYPE=RELEASE .. > cmake.log 
RUN make -j4 > make.log
RUN make install > make_install.log

WORKDIR /

SHELL ["/bin/bash", "-c"]
