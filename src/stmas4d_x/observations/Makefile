# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../../../

# make can't find pgf90, modify makefile.inc and copy to this dir. HJ 4/28/2011
include $(SRCROOT)/src/include/makefile.inc

RM = rm -f

FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a

STINC = -I../

#=============================================================
# STMAS4D obs ingest:
#=============================================================

ARC = libstmas4d_obs.a
DEBUGLIB = libdebugstmas4d_obs.a

SRC= \
    ingest_all.o \
    insitu_height_profiler.o \
    insitu_height_sonde.o \
    insitu_height_surfaces.o \
    back2obs.o \
    quality_control.o
    
ARCHIVE = $(ARC)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all: lib
lib: $(ARC)

$(ARC): $(ARCHIVE)
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib ARC='$(DBLAPSLIBS)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	$(RM) $(EXE)  *.o *~ *# *.mod *.a

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(DBFLAGS) $(FREE) $(STINC) -c $*.f90
.f.o:;      $(FC) $(DBFLAGS) -c $*.f
