# Compiler and compilation flags, using gcc 12 installed by MacPorts
cc=gcc-mp-12
cxx=g++-mp-12 -fopenmp
cflags=-Wall -march=native -ansi -pedantic -O3 -fopenmp -std=c++11 -fPIC -lstdc++fs -Wno-variadic-macros -g -Wuninitialized
