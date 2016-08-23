#!/bin/sh

rm -rf Makefile Makefile.in aclocal.m4 config.status config.h config.h.in config.log configure autom4te.cache config.guess config.sub missing libtool compile depcomp install-sh ltmain.sh

echo "Running aclocal (can take a while)"
aclocal
rc=$?
if [ $rc != 0 ] ; then
  echo
  echo "Oooops, aclocal failed? Maybe a wrong libtool version?"
  exit $rc
fi

echo "Running autoconf"
autoconf
echo "Running autoheader"
autoheader
echo "Running automake"
automake --add-missing

echo "Done. You can now run './configure', then 'make'"
