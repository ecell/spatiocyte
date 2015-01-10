#!/bin/sh
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew doctor
brew prune
brew update
brew tap homebrew/python homebrew/boneyard
brew install Caskroom/cask/xquartz
brew install wget automake autoconf libtool pkg-config gsl pygtk gcc boost-python homebrew/science/hdf5 --with-cxx numpy scipy matplotlib libav glibmm
wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
sudo python ez_setup.py
rm ez_setup.py
sudo rm -rf setuptools*.zip
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
patch -N /usr/local/Cellar/glibmm/2.42.0/include/glibmm-2.4/glibmm.h mac_glibmm_h.diff
brew install gtkglextmm
./autogen.sh
./configure --prefix=$HOME/root
make -j4
make install
