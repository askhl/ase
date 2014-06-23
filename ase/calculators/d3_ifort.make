#!/usr/bin/env sh

f2py --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c d3ef.pyf d3ef.f90
