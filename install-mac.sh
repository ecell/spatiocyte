#!/bin/sh
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew doctor
brew update
brew tap homebrew/python homebrew/boneyard
brew install wget automake autoconf libtool pkg-config gsl pygtk gcc boost-python homebrew/science/hdf5 --with-cxx numpy scipy matplotlib libav gtkglextmm
wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
sudo python ez_setup.py
rm ez_setup.py
wget https://raw.github.com/pypa/pip/master/contrib/get-pip.py
sudo python get-pip.py
rm get-pip.py
sudo pip install ply
echo "alias blender=/Applications/Blender/blender.app/Contents/MacOS/blender" >> ~/.profile
echo "alias vlc=/Applications/VLC.app/Contents/MacOS/VLC" >> ~/.profile
echo "export PATH=$HOME/root/bin:$PATH" >> ~/.profile
echo "export LD_LIBRARY_PATH=$HOME/root/lib:$LD_LIBRARY_PATH" >> ~/.profile
echo "export PYTHONPATH=$HOME/root/lib/python2.7/site-packages:$PYTHONPATH" >> ~/.profile
echo "export PKG_CONFIG_PATH=/opt/X11/lib/pkgconfig" >> ~/.profile
source ~/.profile
cd
mkdir wrk
cd wrk
git clone git://github.com/ecell/spatiocyte
cd spatiocyte
./autogen.sh
./configure --prefix=$HOME/root --disable-gui



vi /usr/local/Cellar/glibmm/2.42.0/include/glibmm-2.4/glibmm.h
patch test.h ori.diff
