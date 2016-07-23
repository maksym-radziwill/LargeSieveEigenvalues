gcc large.c -llapack -lopenblas -o large -Ofast -std=c99
gcc -static large.c -llapack -lopenblas -lgfortran -lquadmath -lm -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -Ofast -o static -std=c99

