#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MSGSIZE 10

int main(int argc, char *argv[]) { 
  MPI_Status stat;
  int numtasks, rank, i, tag = 111, dest = 1, source = 0, count = 0;
  char data[MSGSIZE];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    if (numtasks > 2)
      printf("Number of tasks = %d. Only 2 needed (ignoring extra...)\n", numtasks);
    for (i = 0; i < MSGSIZE; i++)
      data[i] = 'x';
    while (1) {
      MPI_Send(data, MSGSIZE, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
      count ++;
      printf("Sending message number %d to task %d...\n", count, dest);
      sleep(1);
    }
  }
  if (rank == 1) {
    while (1) {
      MPI_Recv(data, MSGSIZE, MPI_CHAR, source, tag, MPI_COMM_WORLD, &stat);
      count++;
      printf("Received message number %d from task %d...\n", count, source);
      sleep(0);
    }
  }
  MPI_Finalize();
}

