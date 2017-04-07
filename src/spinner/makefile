all:  spinner_pg 

spinner_pg : spinner_pg.o readlist_pg.o scaffoldlist_pg.o refine_pg.o gapest_pg.o
	gcc -pg -lm spinner_pg.o readlist_pg.o scaffoldlist_pg.o refine_pg.o gapest_pg.o -o spinner_pg

gapest_pg.o : gapest.c gapest.h scaffoldlist.h readlist.h defs.h
	gcc -pg -c -lm gapest.c -o gapest_pg.o

spinner_pg.o : spinner.c scaffoldlist.h readlist.h defs.h
	gcc -pg -c spinner.c -o spinner_pg.o

readlist_pg.o : readlist.c readlist.h defs.h
	gcc -pg -c readlist.c -o readlist_pg.o

scaffoldlist_pg.o : scaffoldlist.c scaffoldlist.h readlist.h defs.h
	gcc -pg -c scaffoldlist.c -o scaffoldlist_pg.o

refine_pg.o : refine.c refine.h scaffoldlist.h readlist.h defs.h
	gcc -pg -c refine.c -o refine_pg.o

