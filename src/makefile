# Makefile for smis_scaffold 
CC=c++


CFLAGS= -O2  -I${MYSMISDIR}/src/bamtools/include -Wl,-rpath,${MYSMISDIR}/src/bamtools/lib
LFLAGS= -lm -L${MYSMISDIR}/src/bamtools/lib -lbamtools -lz


SMISSORTT= smis_sort.o
SMISRENAME= fast.o smis_rename.o 
SMISSHREDD= fast.o smis_shred.o 
SMISSELCHR= smis_select_chr.o


all  : smis_sort smis_rename smis_shred  smis_select_chr 


smis_sort: makefile $(SMISSORTT)
	$(CC) $(CFLAGS) -o $@ $(SMISSORTT) $(LFLAGS) 
	chmod o-r smis_sort 

smis_select_chr: makefile $(SMISSELCHR)
	$(CC) $(CFLAGS) -o $@ $(SMISSELCHR) $(LFLAGS)
	chmod o-r smis_select_chr

smis_rename: makefile $(SMISRENAME)
	$(CC) $(CFLAGS) -o $@ $(SMISRENAME) $(LFLAGS) 
	chmod o-r smis_rename 

smis_shred: makefile $(SMISSHREDD)
	$(CC) $(CFLAGS) -o $@ $(SMISSHREDD) $(LFLAGS) 
	chmod o-r smis_shred 

