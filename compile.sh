gcc large.c -llapack -lopenblas -lquadmath -o large -march=native -mfpmath=sse -std=c99 -Ofast 
gcc -static large.c -llapack -lopenblas -lgfortran -lquadmath -lm -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -Ofast -o static -std=c99
gcc large_new.c -llapack -lopenblas -lquadmath -o large_float -march=native -mfpmath=sse -std=c99 -Ofast 
gcc -static large_new.c -llapack -lopenblas -lgfortran -lquadmath -lm -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -Ofast -o static_float -std=c99


#-mfpmath=sse -fomit-frame-pointer -floop-parallelize-all -ftree-loop-distribution -ftree-parallelize-loops=4 -fstrict-aliasing -fwhole-program -ffast-math -freciprocal-math 

