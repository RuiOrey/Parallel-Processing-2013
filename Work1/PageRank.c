
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <mpi.h>
#include <stdbool.h>
#include <math.h>

static double *transition_matrix = NULL;
//static double *damping_matrix = NULL;
static double *pagerank_vector = NULL;
static int size_graph = 0;
static int *values_occurrences=NULL;


// since we dynamically allocate our matrix data in a contiguous region
// of memory, we need to use this macro to access a position of the matrix
#define POS(ROW, COL) ((ROW) * size_graph + (COL))



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
use_dummy_graph(double damping)
{
   size_graph = 4;

   const double initial_pagerank = 1.0 / (double)size_graph;
   int i,j;

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

   
   //damping variables
   const double damping_value=1*initial_pagerank*damping;
   const double damping_inverse=1-damping;
   // Generates matrix for dumping - maybe not needed full matrix
   //damping_matrix = (double*)malloc(sizeof(double)*size_graph*size_graph);
   for (i=0;i<size_graph;i++){
      for (j=0;j<size_graph;j++){
         printf("Position %d %d : %f\n",i,j,transition_matrix[POS(i,j)]);
         transition_matrix[POS(i,j)]=transition_matrix[POS(i,j)]*damping_inverse + damping_value; 
         printf("Position %d %d : %f\n",i,j,transition_matrix[POS(i,j)]);
      }
   }

   




}
/*
void generate_transition_matrix(char* line, double damping){
   int i,j, len;
   int num;

   const double initial_pagerank = 1.0 / (double)size_graph;
   const double damping_value=1*initial_pagerank*damping;
   const double damping_inverse=1-damping;

   for (i = 0;i<size_graph;i++)
      for ( j=0;j<size_graph;j++){
         if (values_occurrences[j]!=0)
            transition_matrix[POS(i,j)]= ((1/values_occurrences[j])*damping_inverse)+damping_value;
         else
            transition_matrix[POS(i,j)]= damping_value;

         printf("Position %d %d : %f\n",i,j,transition_matrix[POS(i,j)]);
      }

   }
*/
   void extract_occurrences ( char *line,double damping) {

      int i = 0, len,count=0;
      int line_num,num;
      const double initial_pagerank = 1.0 / (double)size_graph;
      const double damping_value=(double) initial_pagerank*damping;
      const double damping_inverse=(double)1.0-damping;

      sscanf( line, "%d%n", &line_num, &len);
      //printf("Line %d\n",line_num);
      line += len+1;
      while ( sscanf( line, "%d%n", &num, &len) == 1 ) {
         count++;
         transition_matrix[POS(num,line_num)]= 1.0;
         line+=len;
      }
      //printf("line %d : 1/%d\n",line_num,count);
      const float weight= (float) 1.0/count;
      //printf ("wight:%f damping_value:%f",weight,damping_value);
      for(i=0;i<size_graph;i++)
      {
         if(transition_matrix[POS(i,line_num)]==1.0){
            transition_matrix[POS(i,line_num)]=((weight*damping_inverse)+damping_value);
          //  printf("! %f + %f ",(double)weight*damping_inverse,(double)damping_value);
         }
         else{
            transition_matrix[POS(i,line_num)]=(double) (damping_value);
          //  printf("? %f",damping_value);
         }
        // printf(" %f",transition_matrix[POS(i,line_num)]);
      } 
      //printf("\n");
   } 

   static int
   read_file(const char *fname,double damping)
   {
      
      FILE *file;
      if (file = fopen(fname, "r"))
      {    
         int i,j,total;
         //gets size
         char line[1024];
         fgets(line,sizeof(line),file);
         size_graph=atoi(line);

         //calculates page_rank and transition matrix
         const double initial_pagerank = 1.0 / (double)size_graph;
         allocate_pagerank_vector();
         
         for(i = 0; i < size_graph; ++i) {
            pagerank_vector[i] = initial_pagerank;
         }

         printf("%d\n",size_graph);

      // vai ignorar o 1º numero e : e ler para cada linha as ocorrencias dos outros numeros e contar
      //depois vai popular uma matriz em que para cada ocorrencia de um numero mete o valor achado 1/V
      //TODO change line to array of lines

         transition_matrix = (double*)malloc(sizeof(double)*size_graph*size_graph);

         while (fgets(line, sizeof(line), file)) {
            extract_occurrences(line,damping); 
         }



         fclose(file);   
         //sleep(3);
         return 1;

      }
      else{
       printf("File does not exist! Exit\n");

       return 0;
    }
    




   // here we read the file and create the transition matrix and pagerank vector
   // TODO
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

      assert(read_file(file,damping)!=0);
      
      //use_dummy_graph(damping);
      print_transition_matrix();
      print_pagerank_vector();

      // we ensure that the transition matrix can be evenly divided
      assert(size_graph % size == 0);
   } 

   MPI_Bcast(&size_graph, 1, MPI_INT, 0, MPI_COMM_WORLD);
   //MPI_Bcast(&converge, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // non zero processes have not allocated their pagerank vector yet
   if(rank != 0)
      allocate_pagerank_vector();

   // number of rows that each process is getting
   const int rows = size_graph / size;

   // transmit pagerank vector
   int temp_iterator=0;
   double *submatrix = (double *)malloc(sizeof(double)*size_graph * rows);
   double *partial_pagerank_vector = (double *)malloc(sizeof(double) * rows);


   // cycle for iterators 

   int converge=0;
   int converge_temp=0;
   
   while(converge==0)
   { 
      converge_temp=1;
      MPI_Bcast(pagerank_vector, size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /*
      MPI_Scatter involves a designated root process sending data to all processes in a communicator. 
      The primary difference between MPI_Bcast and MPI_Scatter is small but important. 
      MPI_Bcast sends the same piece of data to all processes while MPI_Scatter sends chunks of an 
      array to different processes. Check out the illustration below for further clarification.
      */

      MPI_Scatter(transition_matrix, rows*size_graph, MPI_DOUBLE, submatrix, rows*size_graph,
         MPI_DOUBLE, 0, MPI_COMM_WORLD);
      int numbers=0;
      

      // multiply submatrix
      int i, j;
      for(i = 0; i < rows; ++i) {

         double acc = 0.0;
         for(j = 0; j < size_graph; ++j) {
            acc += submatrix[POS(i, j)] * pagerank_vector[j];
         }
         partial_pagerank_vector[i] = acc;
         converge_temp = converge_temp * (int) (fabs(pagerank_vector[i]-acc)<0.0000009);
         //covergence test - working by reduce
         
         
         // debug messages
         printf("acc:%f item:%d rank:%d before:%f difference:%f boolean :%d converge_temp: %d \n",acc,numbers+rows*rank,rank,pagerank_vector[numbers+rows*rank],fabs(pagerank_vector[numbers+rows*rank]-acc),fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009,converge_temp);
         numbers++;
      }

      /*
      MPI_Gather is the inverse of MPI_Scatter. 
      Instead of spreading elements from one process to many processes, MPI_Gather takes elements 
      from many processes and gathers them to one single process. This routine is highly useful to 
      many parallel algorithms, such as parallel sorting and searching.
      */
      //used for converge
      
      
      MPI_Gather(partial_pagerank_vector, rows, MPI_DOUBLE, pagerank_vector, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      if(rank == 0){
         print_pagerank_vector();
         
      }
      temp_iterator++;
      
      MPI_Barrier(MPI_COMM_WORLD);
      //MPI_Reduce(&converge_temp, &converge, 1,MPI_INT, MPI_PROD,0,MPI_COMM_WORLD);

      if (converge_temp==1) {
         printf("Values Converged at iteration %d!\n",temp_iterator);
         //temp_iterator=iterations; 
      }

   }

   printf("here");
   MPI_Finalize();

   // free dynamic memory
   free(submatrix);
   free(partial_pagerank_vector);
   free(pagerank_vector);
   if(rank == 0)
      {free(transition_matrix);
      }

      return EXIT_SUCCESS;
   }
