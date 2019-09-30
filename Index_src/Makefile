CC=			gcc
CXX=		g++
CFLAGS=		-Wall -g  -O2 
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DMAIN_INDEX -D_FILE_OFFSET_BITS=64  #-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
SRC=		utils.c bwt.c bwtio.c  is.c \
			bntseq.c bwtmisc.c\
		bwt_gen.c  QSufSort.c\
		4bit_bntseq.c 4bit_bwt_gen.c hapmap.c   rbwt.c LookUpTable.c\
	        localPattern.c mixRef.c index1.c	
PROG=		index index_debug
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	. 

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)




index:$(OBJS) 
		$(CC) $(CFLAGS) $(DFLAGS) $(SRC) -o $@ $(LIBS)
	
index_debug:$(OBJS) 
		$(CC) -Wall -g $(DFLAGS) $(SRC) -o $@ $(LIBS)


clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a


