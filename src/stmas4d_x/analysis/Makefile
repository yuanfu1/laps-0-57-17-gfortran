# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../../..

# make can't find pgf90, modify makefile.inc and copy to this dir. HJ 4/28/2011
include $(SRCROOT)/src/include/makefile.inc

RM = rm -f

FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/wind/libwind.a \
           $(SRCROOT)/src/wind/libwind.a  \
           $(SRCROOT)/src/balance/libqbalpe.a \
           $(SRCROOT)/src/lib/temp/libtemp.a \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a \
           $(SRCROOT)/src/stmas4d/normal_subs/libnormal_subs.a \
           $(SRCROOT)/src/stmas4d/observations/libstmas4d_obs.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a \
	     $(SRCROOT)/src/stmas4d/normal_subs/libnormal_subs_debug.a \
	     $(SRCROOT)/src/stmas4d/observations/libstmas4d_obs.a

LAPSFC = $(SRCROOT)/src/sfc/libsfc.a
STINC = -I../

#=============================================================
# STMAS: analysis
#=============================================================

#EXE = STMAS4D.exe
ARC = libstmas4d_anal.a
DEBUGLIB = libdebugstmas4d_anal.a

SRC=STMAS_module.o \
    LAPS_ingest.o \
    get_analysis_window.o \
    grid_interpolation.o \
    STMAS_LAPS_input.o \
    jacobian_cal.o \
    STMAS_analysis.o \
    STMAS_minimizer.o \
    STMAS_costgrad.o \
    STMAS_smooth.o \
    routines.o \
    STMAS_balcost.o \
    mx_costgrad.o \
    my_costgrad.o \
    mz_costgrad.o \
    cont_costgrad.o \
    tv_costgrad.o \
    qv_costgrad.o

ARCHIVE = $(ARC)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all: lib
lib: $(ARC)

$(ARC) : $(ARCHIVE)
	$(RANLIB) $@

debug:  $(SRC)
	$(MAKE) lib ARC='$(DBLAPSLIBS)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)' 

clean:
	$(RM) $(EXE)  *.o *~ *# *.mod

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FFLAGS) $(FREE) $(STINC) -c $*.f90
.f.o:;      $(FC) $(FFLAGS) -c $*.f
