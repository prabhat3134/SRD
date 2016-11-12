# GNU compilers elmo
F95=gfortran
OPTFLAGS=-O2 
#-g -traceback -check all -fp-stack-check  
# FFLAGS=-openmp -I/cvos/shared/apps/openmpi/intel/64/1.2.8/include
FFLAGS = -ffree-form
OBJECTS= SRD_lib.o SRD_main.o   
MODULES= srd_library.mod
EXE=Poiseuille

all: $(EXE)

$(EXE): $(OBJECTS) $(MODULES)
	$(F95) $(FFLAGS) $(OBJECTS) -o $(EXE) $(LINKFLAGS)

%.o:	%.f
	$(F95) $(FFLAGS) -c $<

