
all: carpet2sdf carpet2xgraph

.c.o:
	cc -Wall -g -I$(HOME)/include -c $*.c	

carpet2sdf: carpet2sdf.o
	cc -g carpet2sdf.o -o carpet2sdf -L$(HOME)/lib -lrnpl -lmfhdf -ldf -ljpeg -lz -lvsso -lsv -lm

carpet2xgraph: carpet2xgraph.o
	cc -g carpet2xgraph.o -o carpet2xgraph 
