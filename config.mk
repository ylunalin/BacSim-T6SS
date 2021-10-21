# The compilers used here are GNU GCC and G++ 9 installed with MacPorts.
# Replace the compiler commands and the flags with ones on your computing system.
cc=gcc-mp-9
cxx=g++-mp-9 -fopenmp
cflags=-Wall -march=native -ansi -pedantic -O0 -fopenmp -std=c++11 -fPIC -lstdc++fs -Wno-variadic-macros -g
