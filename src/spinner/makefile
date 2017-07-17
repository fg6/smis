all:  spinner_pg 

spinner_pg : spinner_pg.o readlist_pg.o scaffoldlist_pg.o refine_pg.o gapest_pg.o
	gcc spinner_pg.o readlist_pg.o scaffoldlist_pg.o refine_pg.o gapest_pg.o -o spinner_pg -pg -lm

gapest_pg.o : gapest.c gapest.h scaffoldlist.h readlist.h defs.h 
	gcc -c gapest.c -o gapest_pg.o -pg

spinner_pg.o : spinner.c scaffoldlist.h readlist.h defs.h
	gcc -c spinner.c -o spinner_pg.o -pg

readlist_pg.o : readlist.c readlist.h defs.h
	gcc -c readlist.c -o readlist_pg.o -pg

scaffoldlist_pg.o : scaffoldlist.c scaffoldlist.h readlist.h defs.h
	gcc -c scaffoldlist.c -o scaffoldlist_pg.o -pg

refine_pg.o : refine.c refine.h scaffoldlist.h readlist.h defs.h
	gcc -c refine.c -o refine_pg.o -pg

