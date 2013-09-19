noinst_LTLIBRARIES = libMechanocyte.la
libMechanocyte_la_SOURCES= lupack.f mxlibswD.f iolibswN.f cmlibswD.f interface.f

AM_CPPFLAGS=-I$(srcdir)/../include
AM_CXXFLAGS = -g -Wall -Wextra -O2 \
	-funroll-loops -finline-functions -fomit-frame-pointer

# if HAVE_SSE2
#  AM_CXXFLAGS+= -msse2
# endif
