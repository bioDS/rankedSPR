default: spr_playground.so

spr_playground.so: treeOclock spr_playground.o queue.o
	gcc -shared -o spr_playground.so queue.o ../treeOclock/tree.o ../treeOclock/rnni.o ../treeOclock/spr.o ../treeOclock/unlabelled_rnni.o

spr_playground.o: spr_playground.c spr_playground.h
	gcc -fPIC -Wall -c -g -O2 spr_playground.c

treeOclock:
	make -C ../treeOclock

queue.o: queue.c queue.h
	gcc -fPIC -Wall -c -g -O2 queue.c