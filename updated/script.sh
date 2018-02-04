#!/bin/bash
rm -f *.o
rm -f *.mod
rm -f SRD_serial
rm -f SRD_parallel

if [[ $# -gt 0 ]]
then
	ifort -fpp2  -FR -c -Tf var.f95 
	ifort -fpp2  -FR -c -Tf streaming.f95 
	ifort -fpp2  -FR -c -Tf library_copy.f95
	ifort -fpp2  -FR -o SRD_serial -Tf main_copy.f95 var.o streaming.o library_copy.o
else
	ifort -fpp2  -FR -openmp -c -Tf var.f95 
	ifort -fpp2  -FR -openmp -c -Tf streaming.f95 
	ifort -fpp2  -FR -openmp -c -Tf library_copy.f95
	ifort -fpp2  -FR -openmp -o SRD_parallel -Tf main_copy.f95 var.o streaming.o library_copy.o
fi
