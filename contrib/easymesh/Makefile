include ../../Make.global_options

CC=gcc
LIBS=-lm

default: all
	
all: easymesh

easymesh: easymesh.c
	${CC} -o $@ $^ ${LIBS}

install:
	mkdir -p $(prefix)/bin
	cp easymesh $(prefix)/bin

clean:
	rm -f easymesh

.PHONY: all default clean
