#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  MPI_Status stat;
  int numtasks, rank, tag = 1, alpha, i;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    if (numtasks > 2) 
      printf("Number of tasks = %d. Only 2 needed (ignoring extra...)\n", numtasks);
    for (i = 0; i < 10; i++) {
      alpha = i * 10;
      MPI_Send(&alpha, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
      printf("Task %d sent alpha = %d\n", rank, alpha);
    }
  }
  if (rank == 1) {
    for (i = 0; i < 10; i++) {
      MPI_Recv(&alpha, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
      printf("Task %d received alpha = %d\n", rank, alpha);
    }
  }
  MPI_Finalize();
}
