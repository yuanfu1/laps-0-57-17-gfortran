# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

RM = rm -f
FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/wind/libwind.a \
           $(SRCROOT)/src/lib/temp/libtemp.a \
           $(SRCROOT)/src/lib/liblaps.a  \
	   $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a \
	   $(SRCROOT)/src/sfc/libsfc.a \
	   $(SRCROOT)/src/var/bufr/libbufr.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/wind/libdebugwind.a \
             $(SRCROOT)/src/lib/temp/libdebugtemp.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
	     $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a \
	     $(SRCROOT)/src/sfc/libsfc.a \
             $(SRCROOT)/src/var/bufr/libdebugbufr.a

BUFRLIB = 
DBBUFRLIB = 

HUMIDLIB = $(SRCROOT)/src/humid/liblq3.a

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f

# Prevent make from deleting .o files after building
.PRECIOUS: %.o

.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f
.c.o:;	    $(CC) $(CPPFLAGS) -c $*.c

#=============================================================
# 3DVAR:
#=============================================================

PREP = GSI_Prep.exe
POST = GSI_Post.exe

GSIPR= GSI_Prep.f90 \
       LAPS_Parm.f90 \
       LAPS_Conf.f90 \
       LAPS_Bkgd.f90 \
       LAPS_Obsv.f90 \
       LAPS_Memo.f90 \
       GSI_BkObs.f90 \
       wrfbkgout.f90 \
       StaggerXY.f90 \
       GSI_radarbufr.f90 \
       GSI_noaa1d2bufr.f90 \
       LAPSMassConverter.f90 \
       StaggerConverter.f90 \
       Interpolation.f90

# Object files for GSIPR
GSIPR_OBJ = $(GSIPR:.f90=.o)

GSIPO= GSI_Post.f90 \
       StaggerXY.f90 \
       LAPSMassConverter.f90 \
       StaggerConverter.f90 \
       Interpolation.f90

# Object files for GSIPO
GSIPO_OBJ = $(GSIPO:.f90=.o)

GSICC= getiofile.c

# Object files for GSICC
GSICC_OBJ = $(GSICC:.c=.o)

all: $(PREP) $(POST)

$(PREP): $(GSIPR_OBJ) $(GSICC_OBJ) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(GSIPR_OBJ) $(GSICC_OBJ) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@

$(POST): $(GSIPO_OBJ) $(LAPSLIBS) $(HUMIDLIB)
	$(FC) $(FFLAGS) $(GSIPO_OBJ) $(LAPSFC) $(LAPSLIBS) $(HUMIDLIB) \
	$(BUFRLIB) $(OTHERLIBS) -o $@

debug: $(GSIPR) $(GSIPO)
	$(RM) $(PREP) $(POST);
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS) $(DBBUFRLIB)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: $(PREP) $(POST)
	$(INSTALL) $(PREP) $(INSTALLROOT)/bin/$(PREP)
	$(INSTALL) $(POST) $(INSTALLROOT)/bin/$(POST)

clean:
	rm -f $(PREP) $(POST)  *.o *~ *# *.mod

distclean: clean

