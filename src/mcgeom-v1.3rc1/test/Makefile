# Note that this version DOES NOT USE the OPTIMIZED ATTRIBUTE LIB
#
#
DEVROOT := $(HOME)/develop
ARCHOS := i686_linux

ifeq ($(ARCHOS),i686_linux)
  CC = gcc
  CXX = g++
  CFLAGS = -fPIC -ansi -Wall 
  CXXFLAGS = -fPIC -ansi -Wall 
  CDEPFLAGS = -M
  CXXDEPFLAGS = -M
  RANLIB = ls
  DEBUGFLAGS = -g -DDEBUG
  OPTFLAGS = -O3
  DEFINES = -DLINUX
endif

MCGEOM = ../ 

MCGEOMINC := -I$(MCGEOM)/include

MCGEOMLIB-D = -L$(MCGEOM)/lib/$(ARCHOS) -lmcgeom-d
MCGEOMLIB = -L$(MCGEOM)/lib/$(ARCHOS) -lmcgeom


INCDIR := -I./include $(MCGEOMINC) 

LIBS-D = $(MCGEOMLIB-D) -lm
LIBS = $(MCGEOMLIB) -lm


dirs := src
dirs2 := src
temp := $(foreach dir,$(dirs),$(wildcard $(dir)/*.c))
hdrs := $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
files := $(temp) $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
srcs := $(notdir $(temp))
objs := $(srcs:.c=.o)
obj := $(addprefix obj/$(ARCHOS)-d/,$(objs))
optobj := $(addprefix obj/$(ARCHOS)/,$(objs))
VPATH = obj/$(ARCHOS)-d:obj/$(ARCHOS):$(dirs2)


BINDIR := ./bin/$(ARCHOS)

debug: $(obj)
	$(CC) $(LDFLAGS) $(DEFINES) $(INCDIR) $(DEBUGFLAGS) $^ $(LDOBJS) $(LIBS-D) -o $(BINDIR)/mcgtest-d

obj/$(ARCHOS)-d/%.o: src/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c src/$< -o obj/$(ARCHOS)-d/$@



opt: $(optobj)
	$(CC) $(LDFLAGS) $(DEFINES) $(INCDIR) $(OPTFLAGS) $^ $(LDOBJS) $(LIBS) -o $(BINDIR)/mcgtest

obj/$(ARCHOS)/%.o: src/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c src/$< -o obj/$(ARCHOS)/$@

obj/$(ARCHOS)/%.o: src/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


clean:
	/bin/rm bin/$(ARCHOS)/* obj/$(ARCHOS)-d/*.o obj/$(ARCHOS)/*.o

depend:
	makedepend  $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.
