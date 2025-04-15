#!/bin/bash

clear

# get the pc architecture
MACHINE_TYPE=`uname -m`

key="$1"
# 64-bit compilation
g++ -O3 -std=c++11 -fPIC -c epslib.cpp -lm -Wall  -static-libstdc++
g++ -O3 -std=c++11 -fPIC -c Faddeeva.cc -lm -Wall  -static-libgcc -static-libstdc++
g++ -O3 -shared -o ../clib/epslib-mac.so epslib.o Faddeeva.o -lm -Wall  -static-libstdc++

rm epslib.o
rm Faddeeva.o
