# Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville, Maryland, U.S.A.  All rights reserved.

THISCODEDIR := .
OC1DIR := oc1

# Directories to search for header files
#SEARCHDIRS := -I- -I${THISCODEDIR} -I${OC1DIR}

# C compiler

C       = gcc
CC      = g++
CFLAGS  = -O1 ${SEARCHDIRS}
GFLAGS = -O3 -g -Wall
LIBS    = -lm

OBJ     = util.o tree_util.o load_data.o perturb.o compute_impurity.o impurity_measures.o classify.o prune.o mktree.o
INC     = ${OC1DIR}/oc1.h
MAKEFILE= makefile

%.o : ${OC1DIR}/%.c
	${C} ${CFLAGS} -c ${OC1DIR}/$*.c

.PHONY : all
all:    mktree build1 falsecomp findsites karlin score scoreATG scoreSTOP updomeasure build-imm build-icm dtscore

mktree: $(OBJ)
	${C} $(CFLAGS) $(OBJ) -o mktree $(LIBS)

build1: build1.o
	${CC} build1.c -o build1 $(LIBS)

falsecomp: falsecomp.o
	${CC} falsecomp.c -o falsecomp $(LIBS)

findsites: findsites.o
	${CC} findsites.c -o findsites $(LIBS)

karlin: karlin.o
	${C} karlin.c -o karlin 

score: score.o
	${CC} score.c -o score $(LIBS)

scoreATG: scoreATG.o
	${CC} scoreATG.c -o scoreATG $(LIBS)

scoreSTOP: scoreSTOP.o
	${CC} scoreSTOP.c -o scoreSTOP $(LIBS)

updomeasure: updomeasure.o util2.o
	$(CC) -o updomeasure updomeasure.o util2.o  -lm

updomeasure.o: updomeasure.c 
	$(CC) -c updomeasure.c

util2.o: util2.c
	$(CC) -c util2.c

build-imm: build-imm.o util2.o
	${CC} -o build-imm build-imm.o util2.o $(LIBS)

build-icm: build-icm.o misc.o
	$(C) -o $@ build-icm.o misc.o $(LIBS)

build-icm.o: build-icm.c
	$(C) $(GFLAGS) -c build-icm.c

misc.o: misc.c
	$(C) $(GFLAGS) -c misc.c

dtscore: dtscore.o tree_util_prob.o util2.o
	$(CC) -o dtscore dtscore.o tree_util_prob.o util2.o -lm

dtscore.o: dtscore.c
	$(CC) -c dtscore.c

tree_util_prob.o: tree_util_prob.c
	$(CC) -c tree_util_prob.c
