#!/bin/bash

rm -r *.o

ifort -fpp2 -FR -c -Tf SRD_lib.f95 
ifort -fpp2 -FR -o SRD2D -Tf SRD_main.f95 SRD_lib.o

rm -r *.o


