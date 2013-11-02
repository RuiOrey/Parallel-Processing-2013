#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  MPI_Status stat;
  int numtasks, rank, dest, tag, source, rc, count;
  char inmsg, outmsg='x';

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Task %d starting...\n", rank);
  if (rank == 0) {
    if (numtasks > 2) 
      printf("Number of tasks = %d. Only 2 needed (ignoring extra...)\n", numtasks);
    dest = rank + 1;
    source = rank;
    tag = rank;
    // MPI_Send(&alpha, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
    rc = MPI_Bsend(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    printf("Sent to task %d...\n", dest);
    rc = MPI_Recv(&inmsg, 1, MPI_CHAR, dest, tag+1, MPI_COMM_WORLD, &stat);
    printf("Received from task %d...\n", source);
  } else if (rank == 1) {
    dest = rank - 1;
    source = rank;
    tag = rank;
    // MPI_Recv(&alpha, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);
    rc = MPI_Recv(&inmsg, 1, MPI_CHAR, dest, tag-1, MPI_COMM_WORLD, &stat);
    printf("Received aaaa from task %d...\n", source);
    rc = MPI_Bsend(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    printf("Sent to task %d...\n", dest);
  }
  if (rank < 2) {
    rc = MPI_Get_count(&stat, MPI_CHAR, &count);
    printf("Task %d: received %d char(s) from task %d with tag %d\n",
	   rank, count, stat.MPI_SOURCE, stat.MPI_TAG);
  }
  MPI_Finalize();
}
