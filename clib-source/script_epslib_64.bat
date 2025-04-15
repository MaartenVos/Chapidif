

	g++ -O3 -std=c++11 -fPIC -c  epslib.cpp -lm -Wall -static-libgcc -static-libstdc++ -static 
    g++  -O3 -std=c++11 -fPIC -c Faddeeva.cc -lm -Wall  -static-libgcc -static-libstdc++
	g++ -O3 -shared -o ../clib/epslib-win64.dll  epslib.o Faddeeva.o -lm   -lquadmath -Wall -static-libgcc -static-libstdc++ -static 

del epslib.o
del Faddeeva.o
pause
