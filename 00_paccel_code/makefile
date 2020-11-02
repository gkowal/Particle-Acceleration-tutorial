#-------------------------------------------------------------------------------
#
# Makefile for the GAccel code
#
#-------------------------------------------------------------------------------

HOST    := $(shell hostname -s)
OS      := $(shell uname)

#-------------------------------------------------------------------------------
#
# Default compilation time options
#
#-------------------------------------------------------------------------------

DEBUG         = N
HDF5          = N
RELAT         = N
INTERPOLATION = N
TVD           = N
BNDRY         = N
PERIODIC      = N
PREC          = 8
DIMS          = 3
TEST          = ISLAND
LOSS          = N

# include compilation time options
#
ifeq ($(wildcard make.config),make.config)
$(info Using compilation time options from "$(wildcard make.config)".)
include make.config
else
$(info Using default compilation time options.)
endif

#-------------------------------------------------------------------------------
#
# Default host compiler options
#
#-------------------------------------------------------------------------------

FC            = gfortran
FFLAGS        = -O2
LD            = $(FC)
LDFLAGS       = $(FFLAGS)
LIBS          =

# include the host compilation options
#
ifeq ($(wildcard make.$(HOST)),make.$(HOST))
$(info Using host compiler options from "$(wildcard make.$(HOST))".)
include make.$(HOST)
else
$(info Using default host compiler options.)
endif

#-------------------------------------------------------------------------------
#
# Parse compilation options
#
#-------------------------------------------------------------------------------

ifeq ($(DEBUG),Y)
FFLAGS += -g -DDEBUG
endif
ifeq ($(HDF5),Y)
FFLAGS += -DHDF5
endif
ifeq ($(RELAT),Y)
FFLAGS += -DRELAT
endif
ifeq ($(INTERPOLATION),L)
FFLAGS += -DTRILIN
endif
ifeq ($(INTERPOLATION),C)
FFLAGS += -DTRICUB
ifeq ($(TVD),Y)
FFLAGS += -DTVD
endif
endif
ifeq ($(BNDRY),Y)
FFLAGS += -DBNDRY
endif
ifeq ($(PERIODIC),Y)
FFLAGS += -DPERIODIC
endif
FFLAGS += -DPREC=$(PREC)
FFLAGS += -DDIMS=$(DIMS)
ifeq ($(TEST),WAVE)
FFLAGS += -DTEST -DWTEST
endif
ifeq ($(TEST),ISLAND)
FFLAGS += -DTEST -DITEST
endif
ifeq ($(LOSS),Y)
FFLAGS += -DLOSS
endif

ifeq ($(HDF5),Y)
ifndef $(HDF5DIR)
HDF5DIR = /usr
endif
endif

ifndef $(FITSDIR)
FITSDIR = /usr
endif

ifeq ($(TEST),N)
ifeq ($(HDF5),Y)
# HDF5 include and library paths
#
FFLAGS += -I$(HDF5DIR)/include
LIBS   += -L$(HDF5DIR)/lib
LIBS   += -lhdf5_fortran -lhdf5
endif

# FITS library path
#
LIBS   += -L$(FITSDIR)/lib
LIBS   += -L./ -lcfitsio -lz
endif

#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c $(FFLAGS) $*.F90

#-------------------------------------------------------------------------------

name = paccel

default: $(name).x

sources = fields.F90 fitsio.F90 hdf5io.F90 paccel.F90 params.F90 particles.F90
objects = fields.o   fitsio.o   hdf5io.o   paccel.o   params.o   particles.o
files   = $(sources) makefile params.in license.txt

$(name).x: $(objects)
	$(LD) $(LDFLAGS) $(objects) -o $(name).x $(LIBS)

dist: $(files)
	tar czvf $(name).tar.gz $(files)

clean-bak:
	rm -rf *.bak *~

clean-data:
	rm -rf *.dat info.txt

clean-exec:
	rm -rf *.x

clean-modules:
	rm -rf *.mod

clean-objects:
	rm -rf *.o

clean: clean-bak clean-modules clean-objects

clean-all: clean-bak clean-data clean-exec clean-modules clean-objects

#-------------------------------------------------------------------------------

fields.o        : makefile fields.F90 params.o fitsio.o hdf5io.o
fitsio.o        : makefile fitsio.F90 params.o
hdf5io.o        : makefile hdf5io.F90 params.o
paccel.o        : makefile paccel.F90 fields.o params.o particles.o
params.o        : makefile params.F90
particles.o     : makefile particles.F90 fields.o params.o

#-------------------------------------------------------------------------------
