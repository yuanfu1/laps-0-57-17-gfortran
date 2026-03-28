FC = gfortran
FFLAGS = -O3 -DFORTRANUNDERSCORE
FFLAGS_DB = -g -traceback $(FFLAGS)
FINC = -I./
FLINK=

.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FFLAGS) -c $*.f90
.f.o:;      $(FC) $(FFLAGS) -c $*.f
.c.o:;	    $(CC) $(CFLAGS) -c $*.c

#=============================================================
# Make LAPS mini library:
#=============================================================
REQUIRED= \
          lapack_subs_required.f90 \
          normal_subs.f90

OBJECTS=  \
          lapack_subs_required.o \
          normal_subs.o

#=============================================================
# Build:
#=============================================================
all:
	@echo '';echo '';echo 'Building normal library...';echo '------------------------------'
	$(FC) $(FFLAGS) -c $(REQUIRED) $(FINC) $(FLINK)
	ar r libnormal_subs.a $(OBJECTS)
	ranlib libnormal_subs.a
	@echo '------------------------------';echo '';echo '';echo''

test:
	@echo '';echo '';echo 'Building normal library and test.exe...';echo '------------------------------'
	$(FC) $(FFLAGS) -c $(REQUIRED) $(FINC) $(FLINK)
	ar r libnormal_subs.a $(OBJECTS)
	ranlib libnormal_subs.a
	$(FC) $(FFLAGS) -o TESTNORM.exe testnorm.f90 $(REQUIRED) $(FINC) $(FLINK)
	@echo '------------------------------';echo '';echo '';echo''

debug:
	@echo '';echo '';echo 'Building normal library for debugging...';echo '------------------------------'
	$(FC) $(FFLAGS_DB) -c $(REQUIRED) $(FINC) $(FLINK)
	ar r libnormal_subs_debug.a $(OBJECTS)
	ranlib libnormal_subs_debug.a
	@echo '------------------------------';echo '';echo '';echo''

clean:
	rm -f ./lib/*
	rm -f ./object/*
	rm -f ./modules/*
	rm -f *.o
	rm -f *.mod
	rm -f *~
