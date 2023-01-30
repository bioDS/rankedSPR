default: spr_path.so

spr_path.so: treeOclock seidel spr_path.o
	gcc -shared -o spr_path.so treeOclock/tree.o treeOclock/rnni.o treeOclock/spr.o seidel_compressed/libseidel.so spr_path.o

spr_path.o: spr_path.c spr_path.h
	gcc -fPIC -Wall -c -g -O2 spr_path.c

treeOclock:
	make -C treeOclock

seidel:
	make -C seidel_compressed
