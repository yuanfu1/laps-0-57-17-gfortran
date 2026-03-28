# Used by configure to create makefile.inc 
# 
#
INSTALLROOT=/Users/xiey/developments/da/laps/lapsinstallroot
DATAROOT=/Users/xiey/developments/da/laps/lapsdataroot
CC=gcc

CFLAGS= -DMAC   -Wno-error=int-to-pointer-cast -Wno-error=pointer-to-int-cast -Wno-error=implicit-function-declaration
DBCFLAGS= -DMAC -g -DDEBUG=1



FC       = gfortran
#
# NCAR graphics is used by lapsplot only.
#
NCARGFC  = ncargf90

FFLAGS =  -ffree-line-length-none -fallow-argument-mismatch $(INC) 
LITTLEFLAGS =  $(INC) 
FREE = 
FIXED = 
DBFLAGS =  -ffree-line-length-none $(INC) -g 

# Define pre-processor command line options for .F files.

CPP = gfortran 
CPPFLAGS = $(INC) -I/opt/homebrew/include -Dunknown -DFORTRANUNDERSCORE -DDYNAMIC -DMULTI_EXEC

# Define include paths.
NETCDF=/opt/homebrew
INC = -I$(SRCROOT)/src/include -I$(NETCDF)/include -I$(SRCROOT)/src/lib -I$(SRCROOT)/src/lib/modules -I$(SRCROOT)/src/lib/cloud -I$(SRCROOT)/src/lib/radar/synp -I$(SRCROOT)/src/lib/degrib

AR  = ar
RANLIB   = ranlib
INSTALL=$(SRCROOT)/util/install-sh
PERL = /usr/bin/perl
LIBPATHFLAG=-L

# Grib1 and Grib2 decoder libs and flags
DEGRIBLIBS=-ljasper -lpng -lz 
DEGRIBFLAGS=-DUSE_PNG
CDEGRIBFLAGS=-DMAC -DFORTRANUNDERSCORE

LAPSLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib -llaps
DBLAPSLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib -ldebuglaps

BGDATALIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/bgdata -lbgdata
DBBGDATALIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/bgdata -ldebugbgdata

GRIBLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib -lgrib
DBGRIBLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib -ldebuggrib

MTHERMOLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/mthermo -lmthermo
DBMTHERMOLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/mthermo -ldebugmthermo

UTILLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/util -lutil
DBUTILLIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/util -ldebugutil

OTHERLIBS = $(LIBPATHFLAG)$(NETCDF)/lib -lnetcdff -lnetcdf   -lgcc

GRIB2LIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib2 -lgrib2
DBGRIB2LIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib2 -ldebuggrib2
