SHELL = /bin/sh

INTERNAL_LIB_PATH = ../../lib

INSTALL_PATH = .

GLOBAL_LIBRARIES = \



# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cpp .o .a
.cpp.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 



#LIBS = \
#-L$(INTERNAL_LIB_PATH)/global_utils -lglobal_utils \
#-lm


LIBS = -lm



#LIB_FILES = \
#$(INTERNAL_LIB_PATH)/global_utils/libglobal_utils.a

LIB_FILES = 


#INCLUDES = \
#-I$(INTERNAL_LIB_PATH)/global_utils \
#-I$(INTERNAL_LIB_PATH)/vect_nd


INCLUDES = 

OBJECT_FILES = adjacency_ubin.o bins.o io.o


OBJECT_SRC = adjacency_ubin.cpp bins.cpp io.cpp


RCS_FILES = $(OBJECT_SRC)



#Default target: the offlattice binary
#                "offlattice" is the lattice polymer simulator
adjacency_ubin: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o adjacency_ubin \
	$(OBJECT_FILES) \
	$(LIBS)


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" adjacency_ubin
	cp -f adjacency_ubin $(INSTALL_PATH)

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a *.pure core adjacency_ubin $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/offlattice

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_CPP} -M ${CFLAGS} ${INCLUDES} ${SRCS} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} ${INCLUDES} ${SRCS} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
