#!/usr/bin/env sh

f2py --fcompiler=intelem --f90flags='-openmp' -lmkl_rt -lpthread -lm -liomp5 -c d3ef.pyf d3ef.f90 
