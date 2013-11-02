#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define SIZE 5

int main(int argc, char *argv[]) {
  MPI_Datatype indextype;
  MPI_Status stat;
  int numtasks, rank, errorcode, src = 0, dest, tag = 1, i, blocklengths[2], displacements[2];
  float b[SIZE];
  float a[16] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    if (numtasks != SIZE) {
      printf("Number of tasks must be %d. Aborting...\n", SIZE);
      MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
    
  MPI_Barrier(MPI_COMM_WORLD);

  blocklengths[0] = 2;
  blocklengths[1] = 3;
  displacements[0] = 5;
  displacements[1] = 12;
  
  MPI_Type_indexed(2, blocklengths, displacements, MPI_FLOAT, &indextype);
  MPI_Type_commit(&indextype);

  if (rank == 0)
    for (i = 0; i < numtasks; i++)
      MPI_Send(a, 1, indextype, i, tag, MPI_COMM_WORLD);
  MPI_Recv(b, SIZE, MPI_FLOAT, src, tag, MPI_COMM_WORLD, &stat);
  printf("Task %d vector b = %3.1f %3.1f %3.1f %3.1f %3.1f\n", rank, b[0], b[1], b[2], b[3], b[4]);
  
  MPI_Type_free(&indextype);
  MPI_Finalize();
}
