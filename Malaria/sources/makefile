CC=g++

all: glimmerm

glimmerm: glimmerm.o sites.o
	$(CC) -o glimmerm glimmerm.o sites.o -lm

glimmerm.o: glimmerm.cpp
	$(CC) -c glimmerm.cpp

sites.o: sites.c
	$(CC) -c sites.c
