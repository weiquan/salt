CC=			gcc
CXX=		g++
CFLAGS=		-Wall  -Wno-unused-result -Wno-unused-value -Wno-unused-variable -Wno-unused-function -DHAVE_THREAD -o2 #-pg
CXXFLAGS=	$(CFLAGS)

DFLAGS=		-g -Wall -Wno-unused-result -D_FILE_OFFSET_BITS=64  -DDEBUG#-DHAVE_THREAD ## 

SRC=		bwt.c bwtio.c bwtmisc.c rbwt.c\
		utils.c bntseq.c kstring.c\
		editdistance.c LandauVishkin.c ssw.c\
		variant.c metaref.c lookup.c\
		query.c sam.c\
		aln.c alnse.c alnpe.c indexio.c

OBJS=		bwt.o bwtio.o bwtmisc.o rbwt.o\
		utils.o bntseq.o kstring.o\
		editdistance.o LandauVishkin.o ssw.o\
		variant.o mixRef2.o lookup.o\
		query.o sam.o\
		aln.o alnse.o alnpe.o indexio.o 

PROG=		salt salt_debug 
INCLUDES=		
LIBS=		-lm -lz -lpthread 
SUBDIRS=



all:$(PROG)

salt:$(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $@ $(LIBS)

salt_debug:$(SRC)
	$(CC) $(DFLAGS) $(SRC) -o $@ $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a


