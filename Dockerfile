FROM ubuntu:14.04
RUN apt-get update; apt-get install -y git gfortran libgfortran-4.7-dev autoconf automake libtool g++ libgsl0-dev python-numpy python-ply python-gtk2 libboost-python-dev libgtkmm-2.4-dev libgtkglextmm-x11-1.2-dev libhdf5-serial-dev openssh-server valgrind; git clone git://github.com/ecell/spatiocyte; cd /spatiocyte; git checkout ae931c54c99b6028b7441579b69972706cd8b2cd; ./autogen.sh; ./configure; make; make install
ENV LD_LIBRARY_PATH /usr/local/lib

RUN mkdir /var/run/sshd
RUN echo 'root:screencast' |chpasswd
RUN sed -ri 's/PermitRootLogin without-password/PermitRootLogin yes/g'

EXPOSE 22
CMD    ["/usr/sbin/sshd", "-D"]
