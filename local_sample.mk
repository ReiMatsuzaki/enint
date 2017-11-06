NAEWDY_PATH=${HOME}/src/enint
FC=gfortran
FF=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace -g -cpp -ffree-line-length-512

## -- ifort options --
#FF=-O0 -g -traceback -check all -check bounds -check uninit -debug all -fpp # -DDEBUG_FJSON_PARSER # -DDEBUG_FJSON #-DDEBUG_ISTREAM
#FF+=-DGOOD_MACRO
