# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc




.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .a

# Prevent make from deleting .o files after archiving
.PRECIOUS: %.o

LIB      = ../libbgdata.a
DEBUGLIB = ../libdebugbgdata.a



SRC= namelist_read.F90\
		get_acceptable_grapes.F90 \
		module_configure.f90 \
		module_domain.f90 \
		module_grapes_lga.f90 \
		lga_driver_grapes.F90

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

# Add WRF directory to include path so we can use wrf_lga module
INC_WRF = -I../wrf

.f90.o:
	$(FC) -c $(FFLAGS) $(FREE) $(INC) $(INC_WRF) $*.f90

.F90.o:
	$(FC) -c $(FFLAGS) $(FREE) $(INC) $(INC_WRF) $*.F90

all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

cleanup:
	rm -f *.o *.mod

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean
