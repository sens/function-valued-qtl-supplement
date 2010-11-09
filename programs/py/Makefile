# Standard Makefile for Cython

CC = gcc
CPPFLAGS = -I /usr/lib/python2.6/site-packages/numpy/core/include/ -I /usr/include/python2.6/
CFLAGS = -pipe -march=native -msse4.1 -O2 -Wall -std=c99 -fPIC
CXXFLAGS = $(CFLAGS)
LDFLAGS = -Wl,--as-needed -lm -lpython2.6 -lgsl

%.c : %.pyx
	cython $<

%.so : %.o
	$(CC) -shared -o $@ $< $(LDFLAGS)

.PHONY: clean
clean :
	rm -f *.so *.o

