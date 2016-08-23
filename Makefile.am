SUBDIRS = Src


distclean-local:
	-find . -name ".deps" -exec rm -rf {} \;
	-rm -rf OSXstatlibs OSXstatlibs_bach
	-rm -rf Makefile Makefile.in aclocal.m4 config.status config.h config.h.in config.log configure autom4te.cache config.guess config.sub missing libtool compile depcomp install-sh
	-find . -name "Makefile.in" -exec rm -rf {} \;

clean-local:
	-find . -name "*~" -exec rm -rf {} \;
