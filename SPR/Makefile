default: spr_playground.so

spr_playground.so: treeOclock seidel spr_playground.o queue.o
	gcc -shared -o spr_playground.so queue.o ../treeOclock/tree.o ../treeOclock/rnni.o ../treeOclock/exploring_rnni.o ../treeOclock/spr.o ../seidel/libseidel.so spr_playground.o

spr_playground.o: spr_playground.c spr_playground.h
	gcc -fPIC -Wall -c -g -O2 spr_playground.c

treeOclock:
	make -C ../treeOclock

seidel:
	make -C ../seidel

queue.o: queue.c queue.h
	gcc -fPIC -Wall -c -g -O2 queue.c