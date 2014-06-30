
ARCHOS := i686_linux

ifeq ($(ARCHOS),i686_linux)
  CC = gcc
  CFLAGS = -Wall -ansi -Wno-unused-variable
  DEBUGFLAGS = -g -DDEBUG
  OPTFLAGS = -O
  DEFINES = -DLINUX 
  LDFLAGS := -lm
endif

LIBDIR = ./lib

INCDIR = -I./include 

dirs := src
temp := $(foreach dir,$(dirs),$(wildcard $(dir)/*.c))
hdrs := $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
files := $(temp) $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
srcs := $(notdir $(temp))
objs := $(srcs:.c=.o)
obj-d := $(addprefix obj/$(ARCHOS)-d/,$(objs))
obj := $(addprefix obj/$(ARCHOS)/,$(objs))
VPATH = obj/$(ARCHOS)-d:obj/$(ARCHOS):$(dirs)


# DEBUG VERSION OF THE LIBRARY

debug : $(obj-d)
	 ar rcv $(LIBDIR)/$(ARCHOS)/libmcgeom-d.a $^

obj/$(ARCHOS)-d/%.o: src/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c src/$< -o obj/$(ARCHOS)-d/$@


# OPTIMIZED VERSION OF THE LIBRARY

opt : $(obj)
	ar rcv $(LIBDIR)/$(ARCHOS)/libmcgeom.a $^ 

obj/$(ARCHOS)/%.o: src/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c src/$< -o obj/$(ARCHOS)/$@


clean :
	rm $(obj-d) $(obj) $(LIBDIR)/$(ARCHOS)/*.a

print: $(files)
	ppc $?
	touch print

printh: $(hdrs)
	ppc $?
	touch  printh

depend:
	makedepend $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.
