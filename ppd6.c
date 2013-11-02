#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define REQS 1000
#define DISP 100

int main(int argc, char *argv[]) {
  int numtasks, rank, buf, tag = 1, i, errorcode, src, dest, offset, nreqs;
  MPI_Request reqs[REQS * 4];
  MPI_Status stats[REQS * 2];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    if (numtasks != 4) {
      printf("Number of tasks must be 4. Aborting...\n");
      MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
<<<<<<< HEAD
  MPI_Barrier(MPI_COMM_WORLD); 
=======
  MPI_Barrier(MPI_COMM_WORLD);
>>>>>>> d6a01d89615f21ed1251eeacc059a281a73ff75d

  /* tasks 0 and 1 do a Isend/Irecv test */
  if (rank < 2) {
    nreqs = REQS * 2;
    if (rank == 0) {
      src = 1;
      offset = 0;
    }
    if (rank == 1) {
      src = 0;
      offset = REQS * 2;
    }
    dest = src;
    for (i = 0; i < REQS; i++) {
      MPI_Isend(&rank, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &reqs[offset]);
<<<<<<< HEAD
      MPI_Irecv(&buf, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &reqs[offset ]);
      offset += 2;
      if ((i+1) % DISP == 0)
	printf("Task %d has done %d isends/irecvs\n", rank, i + 1);
    }
  }
  
  /* tasks 2 and 3 do a Send/Irecv test */
  if (rank > 1) {
  	
=======
      MPI_Irecv(&buf, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &reqs[offset + 1]);
      offset += 2;
      if ((i+1) % DISP == 0)
  printf("Task %d has done %d isends/irecvs\n", rank, i + 1);
    }
  }

  /* tasks 2 and 3 do a Send/Irecv test */
  if (rank > 1) {
>>>>>>> d6a01d89615f21ed1251eeacc059a281a73ff75d
    nreqs = REQS;    
    if (rank == 2) {
      dest = 3;
      for (i = 0; i < REQS; i++) {
<<<<<<< HEAD
	MPI_Send(&rank, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		
	if ((i+1) % DISP == 0)
	  printf("Task %d has done %d sends\n", rank, i+1);
=======
  MPI_Send(&rank, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
  MPI_Wait(reqs,stats);
  
  if ((i+1) % DISP == 0)
    printf("Task %d has done %d sends\n", rank, i+1);
>>>>>>> d6a01d89615f21ed1251eeacc059a281a73ff75d
      }
    }
    if (rank == 3) {
      src = 2;
      offset = 0;
<<<<<<< HEAD
      for (i = 0; i < REQS-1; i++) {
		MPI_Irecv(&buf, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &reqs[offset]);
		 MPI_Wait(&reqs[offset],MPI_STATUS_IGNORE);
		offset += 1;

		if ((i+1) % DISP == 0)
		  printf("Task %d has done %d irecvs\n", rank, i+1);
		  }
		} 
   
=======
      for (i = 0; i < REQS; i++) {
  MPI_Irecv(&buf, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &reqs[offset]);
  
  
  offset += 1;
  if ((i+1) % DISP == 0)
    printf("Task %d has done %d irecvs\n", rank, i+1);
      }
    } 
  
>>>>>>> d6a01d89615f21ed1251eeacc059a281a73ff75d
  }

  MPI_Waitall(nreqs, reqs, stats);
  MPI_Finalize();
<<<<<<< HEAD
}
=======
}
>>>>>>> d6a01d89615f21ed1251eeacc059a281a73ff75d
