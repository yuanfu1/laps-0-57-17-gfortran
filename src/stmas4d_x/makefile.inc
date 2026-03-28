# Used by configure to create makefile.inc 
# 
#
INSTALLROOT=/usr/nfs/lapb/parallel/laps
DATAROOT=/data/lapb/parallel/laps/data
CC=cc

CFLAGS= -DLITTLE -DSWAPBYTE  -O2 
DBCFLAGS= -DLITTLE -DSWAPBYTE -g -DDEBUG=1



FC       = /opt/pgi/linux86-64/bin/pgf90
#
# NCAR graphics is used by lapsplot only.
#
NCARGFC  = ncargf90

FFLAGS =  -I$(SRCROOT)/src/lib -Mnofree -byteswapio -i4 $(INC) -O2 -fastsse -Bstatic
LITTLEFLAGS =  -I$(SRCROOT)/src/lib -Mnofree -i4 $(INC) -O2 -fastsse
FREE = -Mfree
FIXED = -Mnofree
DBFLAGS =  -I$(SRCROOT)/src/lib -Mnofree -byteswapio -i4 $(INC) -g -Mbounds 

# Define pre-processor command line options for .F files.

CPP = /opt/pgi/linux86-64/bin/pgf90
CPPFLAGS = $(INC)  -Dx86_64 -DF90 -DFORTRANUNDERSCORE -DDYNAMIC -DMULTI_EXEC

# Define include paths.
NETCDF=/opt/nfs/netcdf
INC = -I$(SRCROOT)/src/include -I$(NETCDF)/include -I$(SRCROOT)/src/lib -I$(SRCROOT)/src/lib/modules -I$(SRCROOT)/src/lib/radar/synp

AR  = ar
RANLIB   = ranlib
INSTALL=$(SRCROOT)/util/install-sh
PERL = /usr/bin/perl
LIBPATHFLAG=-L

# Grib1 and Grib2 decoder libs and flags
DEGRIBLIBS=-ljasper -lpng -lz 
DEGRIBFLAGS=-DUSE_JPEG2000 -DUSE_PNG
CDEGRIBFLAGS=-DLINUX -DSWAPBYTE -DFORTRANUNDERSCORE

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

OTHERLIBS = $(LIBPATHFLAG)$(NETCDF)/lib -lnetcdf  

GRIB2LIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib2 -lgrib2
DBGRIB2LIB=$(LIBPATHFLAG)$(SRCROOT)/src/lib/grib2 -ldebuggrib2
