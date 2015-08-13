FROM ubuntu:14.04
RUN apt-get update; apt-get install -y x11vnc xvfb xfce4 tightvncserver git autoconf automake libtool g++ libgsl0-dev python-numpy python-ply python-gtk2-dev libboost-dev libboost-python-dev libgtkmm-2.4-dev libgtkglextmm-x11-1.2-dev libhdf5-dev openssh-server valgrind; git clone git://github.com/ecell/spatiocyte; cd /spatiocyte; ./autogen.sh; ./configure; make; make install; cd; echo "export LD_LIBRARY_PATH=/usr/local/lib" >> .bashrc
RUN mkdir ~/.vnc; touch ~/.vnc/passwd
RUN x11vnc -storepasswd "ecell" ~/.vnc/passwd
ENV USER root
