CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall #-O2#-pg 
CXXFLAGS=	$(CFLAGS)

DFLAGS=		-DHAVE_PTHREAD#-DDEBUG#-D_FILE_OFFSET_BITS=64 

DFLAGS=		#-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
OBJS=		utils.o lv.o bntseq.o kstring.o samParser.o ssw.o polish.o
PROG=		polish
INCLUDES=		
LIBS=		-lm -lz #-lpthread 
SUBDIRS=	. 

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

polish:$(OBJS) 
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)
		rm -f *.o
clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a


