# Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville, Maryland, U.S.A.  All rights reserved.



CC=g++ -g

all: glimmerm

glimmerm: glicmerm.o sites.o tree_util_prob.o util.o
	$(CC) -o glimmerm glicmerm.o sites.o tree_util_prob.o util.o  -lm 

glicmerm.o: glicmerm.cpp 
	$(CC) -c glicmerm.cpp

sites.o: sites.c
	$(CC) -c sites.c

tree_util_prob.o: tree_util_prob.c
	$(CC) -c tree_util_prob.c

util.o: util.c
	$(CC) -c util.c
