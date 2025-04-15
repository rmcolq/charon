FROM ubuntu:20.04

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    make \
    build-essential \
    cmake \
    git

RUN mkdir -p /charon/build
WORKDIR /charon/build

RUN cmake -D CMAKE_BUILD_TYPE=RELEASE .. &> cmake.log
RUN make -j4 &> make.log
RUN make install &> make-install.log

ENV PATH=/charon/build/bin/:$PATH

CMD ["/bin/bash"]
