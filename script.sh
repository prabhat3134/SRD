#!/bin/bash

rm -f *.o
rm -f *.mod
if [ $# -gt 1 ]; then
	ifort -fpp2  -FR -c -Tf SRD2D_ifort_library.f95
	ifort -fpp2  -FR -o SRD2D -Tf SRD2D_main.f95 SRD2D_ifort_library.o
elif [ $# -gt 0 ]; then
	gfortran -ffree-line-length-none -fopenmp -c SRD2D_gcc_library.f95 -fcheck=bounds
	gfortran -ffree-line-length-none -fopenmp -o SRD2D_par SRD2D_main.f95 SRD2D_gcc_library.o -fcheck=bounds
else
	gfortran -ffree-line-length-none -cpp -c SRD2D_gcc_library.f95 -fcheck=bounds
	gfortran -ffree-line-length-none -cpp -o SRD2D_serial SRD2D_main.f95 SRD2D_gcc_library.o -fcheck=bounds
fi
