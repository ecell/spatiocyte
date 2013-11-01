.NOTPARALLEL:
noinst_LTLIBRARIES = libMechanocyte.la
libMechanocyte_la_SOURCES= lupack.f iolibswN.f mxlibswD.f cmlibswD.f molibswE.f avlibsw45.f deformlib.f interface.f

AM_CPPFLAGS=-I$(srcdir)/../include
AM_CXXFLAGS = -g -Wall -Wextra -O2 \
	-funroll-loops -finline-functions -fomit-frame-pointer
FFLAGS = -O3

# if HAVE_SSE2
#  AM_CXXFLAGS+= -msse2
# endif
