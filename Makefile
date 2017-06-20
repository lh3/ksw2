CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat 
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=	-I.
OBJS=		ksw2_gg.o ksw2_gg2.o ksw2_gg2_sse.o ksw2_gg2_sse_u.o \
			ksw2_extz.o ksw2_extz2_sse.o ksw2_extz2_sse_u.o
PROG=		ksw2-test
LIBS=		-lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ksw2-test:cli.o kalloc.o $(OBJS)
		$(CC) $(CFLAGS) $^ -o $@  $(LIBS)
		
clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

cli.o: ksw2.h kseq.h
kalloc.o: kalloc.h
ksw2_extz.o: ksw2.h
ksw2_extz2_sse.o: ksw2.h
ksw2_extz2_sse_u.o: ksw2.h
ksw2_gg.o: ksw2.h
ksw2_gg2.o: ksw2.h
ksw2_gg2_sse.o: ksw2.h
ksw2_gg2_sse_u.o: ksw2.h
