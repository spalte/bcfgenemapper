PROG=		bcfgenemapper

all: $(PROG)
build: all

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
DFLAGS=
OBJS=		main.o genemapper.o
INCLUDES=	-I. -I$(HTSDIR)

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
mandir      = $(prefix)/share/man
man1dir     = $(mandir)/man1

INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644


all:$(PROG) plugins

# See htslib/Makefile
PACKAGE_VERSION  = 0.0.1
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define BCFGENEMAPPER_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all distclean install lib tags test testclean force plugins

force:

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

main.o: version.h $(HTSDIR)/version.h
genemapper.o: genemapper.c genemapper.h

bcfgenemapper: $(HTSLIB) $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm -ldl


clean:
		rm -fr *.o *.dSYM *~ $(PROG) version.h

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch]
