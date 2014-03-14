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
	echo '#define VCFPARSER_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all distclean install lib tags test testclean force plugins

force:

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

test: $(PROG) plugins test/test-rbuf
		./test/test.pl

PLUGINC = $(foreach dir, plugins, $(wildcard $(dir)/*.c))
PLUGINS = $(PLUGINC:.c=.so)

plugins: $(PLUGINS)

%.so: %.c vcfannotate.c $(HTSDIR)/libhts.so
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ -L$(HTSDIR) $< -lhts

main.o: version.h $(HTSDIR)/version.h
genemapper.o: genemapper.c genemapper.h

bcfgenemapper: $(HTSLIB) $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm -ldl


install: $(PROG)
		mkdir -p $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
		$(INSTALL_PROGRAM) $(PROG) plot-vcfstats vcfutils.pl $(DESTDIR)$(bindir)
		$(INSTALL_DATA) bcftools.1 $(DESTDIR)$(man1dir)

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG) version.h plugins/*.so

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
