default: treeOclock
	# gcc -fPIC -Wall -c -g -O2 -fsanitize=address tree.c
	gcc -fPIC -Wall -c -g -O2 spr_playground.c
	gcc -shared -o spr_playground.so ../treeOclock/tree.o ../treeOclock/rnni.o ../treeOclock/spr.o ../treeOclock/queue.o ../treeOclock/unlabelled_rnni.o

treeOclock:
	make -C ../treeOclock