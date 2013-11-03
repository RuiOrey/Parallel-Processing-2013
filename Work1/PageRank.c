
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <mpi.h>

static double *transition_matrix = NULL;
static double *pagerank_vector = NULL;
static int size_graph = 0;

// since we dynamically allocate our matrix data in a contiguous region
// of memory, we need to use this macro to access a position of the matrix
#define POS(ROW, COL) ((ROW) * size_graph + (COL))

static void
read_file(const char *file)
{
   // here we read the file and create the transition matrix and pagerank vector
   // TODO
}

static void
print_transition_matrix(void)
{
   int i, j;

   assert(size_graph > 0);

   printf("\nTransition Matrix:\n");
   for(i = 0; i < size_graph; ++i) {
      printf("\t");
      for(j = 0; j < size_graph; ++j) {
         if(j > 0)
            printf(" ");
         printf("%lf", transition_matrix[POS(i, j)]);
      }
      printf("\n");
   }
}

static void
print_pagerank_vector(void)
{
   int i;

   assert(size_graph > 0);
   printf("\nPageRank Vector:\n\t");
   for(i = 0; i < size_graph; ++i) {
      if(i > 0)
         printf(" ");
      printf("%lf", pagerank_vector[i]);
   }
   printf("\n");
}

static void
allocate_pagerank_vector(void)
{
   assert(size_graph > 0);
   assert(pagerank_vector == NULL);
   pagerank_vector = (double*)malloc(sizeof(double)*size_graph);
}

static void
use_dummy_graph(void)
{
   size_graph = 4;

   const double initial_pagerank = 1.0 / (double)size_graph;
   int i;

   allocate_pagerank_vector();
   transition_matrix = (double*)malloc(sizeof(double)*size_graph*size_graph);

   for(i = 0; i < size_graph; ++i) {
      pagerank_vector[i] = initial_pagerank;
   }

   transition_matrix[POS(0, 0)] = 0.0;
   transition_matrix[POS(0, 1)] = 0.0;
   transition_matrix[POS(0, 2)] = 1.0;
   transition_matrix[POS(0, 3)] = 1.0/2.0;
   transition_matrix[POS(1, 0)] = 1.0/3.0;
   transition_matrix[POS(1, 1)] = 0.0;
   transition_matrix[POS(1, 2)] = 0.0;
   transition_matrix[POS(1, 3)] = 0.0;
   transition_matrix[POS(2, 0)] = 1.0/3.0;
   transition_matrix[POS(2, 1)] = 1.0/2.0;
   transition_matrix[POS(2, 2)] = 0.0;
   transition_matrix[POS(2, 3)] = 1.0/2.0;
   transition_matrix[POS(3, 0)] = 1.0/3.0;
   transition_matrix[POS(3, 1)] = 1.0/2.0;
   transition_matrix[POS(3, 2)] = 0.0;
   transition_matrix[POS(3, 3)] = 0.0;
}

int
main(int argc, char **argv)
{
   if(argc != 5) {
      fprintf(stderr, "usage: ./pagerank <graph file> <damping factor> <iterations> <method>\n");
      return EXIT_FAILURE;
   }

   const char *file = argv[1];
   const double damping = atof(argv[2]);
   const int iterations = atoi(argv[3]);
   const char *method = argv[4];

   if(damping <= 0.0 || damping >= 1.0) {
      fprintf(stderr, "Damping factor must be between 0.0 and 1.0\n");
      return EXIT_FAILURE;
   }
   if(iterations < 1) {
      fprintf(stderr, "Number of iterations must be positive\n");
      return EXIT_FAILURE;
   }

   int rank, size;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if(rank == 0) {
      printf("Rank: %d\n", rank);
      printf("Size: %d\n", size);
      printf("File: %s\n", file);
      printf("Damping: %lf\n", damping);
      printf("Iterations: %d\n", iterations);

      // read_file(file);
      use_dummy_graph();
      print_transition_matrix();
      print_pagerank_vector();

      // we ensure that the transition matrix can be evenly divided
      assert(size_graph % size == 0);
   }

   MPI_Bcast(&size_graph, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // non zero processes have not allocated their pagerank vector yet
   if(rank != 0)
      allocate_pagerank_vector();

   // number of rows that each process is getting
   const int rows = size_graph / size;

   // transmit pagerank vector
   MPI_Bcast(pagerank_vector, size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   double *submatrix = (double *)malloc(sizeof(double)*size_graph * rows);
   double *partial_pagerank_vector = (double *)malloc(sizeof(double) * rows);

   MPI_Scatter(transition_matrix, rows*size_graph, MPI_DOUBLE, submatrix, rows*size_graph,
         MPI_DOUBLE, 0, MPI_COMM_WORLD);

   // multiply submatrix
   int i, j;
   for(i = 0; i < rows; ++i) {
      double acc = 0.0;
      for(j = 0; j < size_graph; ++j) {
         acc += submatrix[POS(i, j)] * pagerank_vector[j];
      }
      partial_pagerank_vector[i] = acc;
   }

   MPI_Gather(partial_pagerank_vector, rows, MPI_DOUBLE, pagerank_vector, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   if(rank == 0)
      print_pagerank_vector();

   MPI_Finalize();

   // free dynamic memory
   free(submatrix);
   free(partial_pagerank_vector);
   free(pagerank_vector);
   if(rank == 0)
      free(transition_matrix);

   return EXIT_SUCCESS;
}
