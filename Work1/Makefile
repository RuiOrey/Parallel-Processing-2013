all: compile

compile:
	mpicc -c PageRank.c
	mpicc -o PageRank PageRank.o
	
run: compile
	mpirun -np 4 -hostfile ~/Datasets/dummy.graph 0.15 20 notfox
