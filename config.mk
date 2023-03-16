# Compiler and compilation flags, using gcc 12 installed by MacPorts
cc=gcc
cxx=g++ -fopenmp
cflags=-Wall -march=native -ansi -pedantic -O3 -fopenmp -std=c++11 -fPIC -lstdc++fs -Wno-variadic-macros -g -Wuninitialized
