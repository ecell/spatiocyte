#!/bin/sh -x
sudo apt-get install git autoconf automake libtool g++ libgsl0-dev python-numpy python-ply python-gtk2 libboost-python-dev libgtkmm-2.4-dev libgtkglextmm-x11-1.2-dev libhdf5-serial-dev libav-tools blender vlc python-numpy python-scipy python-matplotlib valgrind
echo "export PATH=.:\$HOME/root/bin:\$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$HOME/root/lib:\$LD_LIBRARY_PATH:." >> ~/.bashrc
echo "export PYTHONPATH=\$HOME/root/lib/python2.7/site-packages:\$PYTHONPATH" >> ~/.bashrc
echo "export ECELL3_DM_PATH=." >> ~/.bashrc
source ~/.bashrc
cd
mkdir wrk
cd wrk
git clone git://github.com/ecell/spatiocyte
cd spatiocyte
./autogen.sh
./configure --prefix=$HOME/root
make -j4
make install
echo "Close and reopen this terminal before proceeding"
