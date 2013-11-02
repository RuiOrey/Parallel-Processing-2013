#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define ARRAYSIZE 160
#define MASTER	  0
float data[ARRAYSIZE];
float do_work(int rank, int offset, int chunksize);

int main(int argc, char *argv[]) {
  int numtasks, rank, errorcode, i, j, tag1, tag2, src, dest, offset, chunksize; 
  float mysum, sum;
  MPI_Status status;
   
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    if (numtasks % 4 != 0) {
      printf("Number of tasks must be divisible by 4. Aborting...\n");
      MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  chunksize = (ARRAYSIZE / numtasks);
  tag1 = 2;
  tag2 = 1;  
  if (rank == MASTER) {
    sum = 0;
    for (i = 0; i < ARRAYSIZE; i++) {
      data[i] = i * 1.0;
      sum = sum + data[i];
    }
    printf("Initial sum = %f\n", sum);
    offset = chunksize;
    for (dest = 1; dest < numtasks; dest++) {
      MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
      MPI_Send(&data[offset], chunksize, MPI_FLOAT, dest, tag2, MPI_COMM_WORLD);
      printf("Sent %d elements to task %d offset = %d\n", chunksize, dest, offset);
      offset = offset + chunksize;
    }
    offset = 0;
    mysum = do_work(rank, offset, chunksize);
    sum = mysum;
    for (i = 1; i < numtasks; i++) {
      src = i;
      MPI_Recv(&offset, 1, MPI_INT, src, tag1, MPI_COMM_WORLD, &status);
      MPI_Recv(&data[offset], chunksize, MPI_FLOAT, src, tag2, MPI_COMM_WORLD, &status);
            //MPI_Recv(&mysum, 1, MPI_FLOAT, src, MPI_SUM, MPI_COMM_WORLD, &status);
      //sum += mysum;
    }
    MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    printf("\nResults\n");
    offset = 0;
    for (i = 0; i < numtasks; i++) {
      for (j = 0; j < chunksize; j++) 
	printf("%f ", data[offset + j]);
      printf("\n");
      offset = offset + chunksize;
    }
    printf("\nFinal sum = %f\n", sum);
  }

  if (rank != MASTER) {
    src = MASTER;
    MPI_Recv(&offset, 1, MPI_INT, src, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&data[offset], chunksize, MPI_FLOAT, src, tag2, MPI_COMM_WORLD, &status);
    mysum = do_work(rank, offset, chunksize);
    dest = MASTER;
    MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
    MPI_Send(&data[offset], chunksize, MPI_FLOAT, dest, tag2, MPI_COMM_WORLD); 
    MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, dest, MPI_COMM_WORLD);
  }

  MPI_Finalize(); 
}


float do_work(int rank, int offset, int chunksize) {
  float sum;
  int i; 

  sum = 0;
  for (i = offset; i < offset + chunksize; i++) {
    data[i] = data[i] * 2;
    sum = sum + data[i];
  }
  printf("Task %d sum = %f\n", rank, sum);
  return sum;
}
