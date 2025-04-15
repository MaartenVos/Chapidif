#!/bin/bash

clear
mkdir -p ../clib
# get the pc architecture
MACHINE_TYPE=`uname -m`

key="$1"

if [ ${MACHINE_TYPE} == 'x86_64' ]; then

    g++  -O3 -std=c++11 -fPIC -c epslib.cpp -Wall  -static-libgcc -static-libstdc++
    g++  -O3 -std=c++11 -fPIC -c Faddeeva.cc   -Wall  -static-libgcc -static-libstdc++
#clang  -Wall -Wextra  -fvisibility=default  -pedantic -c -fPIC epslib.cpp  -o epslib.o 
	g++  -O3 -shared -o ../clib/epslib-lin64.so epslib.o  Faddeeva.o -lm -lquadmath -Wall -static-libgcc #-static-libstdc++
#clang -shared -fvisibility=default -o ../clib/epslib-lin64.so epslib.o  -lm -Wall -static-libgcc
# clang++ -o epslib.os -c -std=c++11 -Wall -fvectorize -fslp-vectorize -fcolor-diagnostics -O2 -fPIC epslib.cpp
# clang++ -o epslib.o -shared epslib.os 
else
	g++  -O3 -std=c++11 -fPIC -c epslib.cpp -lm -Wall  -static-libgcc -static-libstdc++
    g++  -O3 -std=c++11 -fPIC -c Faddeeva.cc -lm -Wall  -static-libgcc -static-libstdc++
	g++   -O3 -shared -o ../clib/epslib-lin32.so epslib.o Faddeeva.o -lm -Wall -static-libgcc -static-libstdc++
fi

rm epslib.o
rm Faddeeva.o
#add -g for debugging
