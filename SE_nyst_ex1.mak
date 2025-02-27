# Makefile ( n source file 1 execute file version )

PACKAGE	= SE_nyst_ex1
SRCS	= $(PACKAGE).c SE_trans.c matrixvector.c sici.c polevl.c
HEADS	= SE_trans.h sigma.h matrixvector.h mconf.h
OBJS	= $(SRCS:.c=.o)
FILES	= $(PACKAGE).mak $(HEADS) $(SRCS)
VER	= `date +%Y%m%d`


### command and flags ###
# uncomment when debugging
#DEBUG	= -ggdb -pg # -lefence

# common (*.o)
LD	= gcc
LDFLAGS	= -g $(DEBUG)
LDLIBS	= -framework Accelerate

# C (*.c)
CC	= gcc
OPTIMIZE= #-O2
CFLAGS	= -g $(OPTIMIZE) -Wall $(DEBUG)
CPPFLAGS= -I.

# C++ (*.cpp)
CXX	= g++
CXXFLAGS= -g $(OPTIMIZE) -Wall $(DEBUG)

# Fortran77 (*.f)
FC	= f77
FFLAGS	= -Wall $(DEBUG)

# Pascal (*.p)
PC	= pc
PFLAGS	= -Wall $(DEBUG)

# etc
SHELL	= /bin/sh
RM	= rm -f
PROF	= gprof


### rules ###

.SUFFIXES:
.SUFFIXES: .o .c .cpp .f .p

all: $(PACKAGE)

$(PACKAGE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $@ $(LDLIBS)

$(OBJS): $(HEADS) $(PACKAGE).mak

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
.f.o:
	$(FC) $(FFLAGS) -c $< -o $@
.p.o:
	$(PC) $(PFLAGS) $(CPPFLAGS) -c $< -o $@


### useful commands ###

clean:
	$(RM) $(PACKAGE) $(OBJS)
	$(RM) core gmon.out *~ #*#

tar:
	@echo $(PACKAGE)-$(VER) > .package
	@$(RM) -r `cat .package`
	@mkdir `cat .package`
	@ln $(FILES) `cat .package`
	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
	@$(RM) -r `cat .package` .package

zip:
	zip -9 $(PACKAGE)-$(VER).zip $(FILES)


prof: run
	$(PROF) $(PACKAGE) | less

run: all
#	./$(PACKAGE) < sample-data | less
	./$(PACKAGE)
