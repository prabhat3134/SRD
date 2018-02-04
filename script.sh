#!/bin/bash

rm -f *.o
rm -f *.mod

if [[ $# -gt 0 ]]
then
	ifort -fpp2  -FR -c -Tf SRD2D_ifort_library.f95
	ifort -fpp2  -FR -o SRD2D -Tf SRD2D_main.f95 SRD2D_ifort_library.o
else
	# gfortran -ffree-line-length-none -c SRD2D_gcc_library.f95 -fcheck=bounds
	# gfortran -ffree-line-length-none -o SRD2D -c SRD2D_main.f95 SRD2D_gcc_library.o -fcheck=bounds
	gfortran -ffree-line-length-none -c SRD2D_gcc_library.f95
	gfortran -ffree-line-length-none -o SRD2D -c SRD2D_main.f95 SRD2D_gcc_library.o 
fi
