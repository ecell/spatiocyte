FROM ubuntu:14.04
RUN apt-get update; apt-get install -y git gfortran libgfortran-4.7-dev autoconf automake libtool g++ libgsl0-dev python-numpy python-ply libboost-python-dev libgtkmm-2.4-dev libgtkglextmm-x11-1.2-dev libhdf5-serial-dev valgrind

