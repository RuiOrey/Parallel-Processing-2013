
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <mpi.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>


static double *transition_matrix = NULL;
//static double *damping_matrix = NULL;
static double *pagerank_vector = NULL;
static int size_graph = 0;
static int *values_occurrences=NULL;
static int converge;
static int converge_temp=1;
static int temp_iterator;
static int method_flag=0;
static int q=0;
static int n_fox=0;
MPI_Comm new_comm,row,col;
MPI_Datatype submatrix_type;

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

void reorder_matrix(int size){
 double *transition_matrix_temp = (double*)malloc(sizeof(double)*size_graph*size_graph);
 int k;
 int line_actual,col_actual,q_actual_linha,k_counter;
 q_actual_linha=0;                  //linha de processos
 k_counter=0;                       //para frequencia de processos
 int q_actual_coluna=0;             //coluna de processos
 int position=0;
 int colu;

 int start_col_index;

 printf("Enter reorder\n");
 for (k=0;k<size;k++)
 {

   if(k_counter>q-1){
      q_actual_linha=q_actual_linha+n_fox;
      k_counter=0;
   }

   q_actual_coluna=abs(q_actual_linha*q-k);
   printf("q%d k_counter%d\n",q,k_counter);

   for (line_actual=0;line_actual<n_fox;line_actual++){
      for (col_actual=0;col_actual<n_fox;col_actual++){

        colu=k_counter*n_fox+col_actual;

        transition_matrix_temp[position] = transition_matrix[POS(line_actual+q_actual_linha,colu)];
        position++;

        printf(" Process:%d and matrix absolute Position:%d will have contents of line %d and column %d of the old matrix\n",k,position,line_actual+q_actual_linha,colu);
     }
  }
  k_counter++;
}
memcpy(transition_matrix,transition_matrix_temp,sizeof(double)*size_graph*size_graph);

}


void reorder_matrix_correctly(int size){
 double *transition_matrix_temp = (double*)malloc(sizeof(double)*size_graph*size_graph);
 int k;
 int line_actual,col_actual,q_actual_linha,k_counter;
 q_actual_linha=0;                  //linha de processos
 k_counter=0;                       //para frequencia de processos
 int q_actual_coluna=0;             //coluna de processos
 int position=0;
 int colu;

 int start_col_index;

 printf("Enter reorder\n");
 for (k=0;k<size;k++)
 {

   if(k_counter>q-1){
      q_actual_linha=q_actual_linha+n_fox;
      k_counter=0;
   }

   q_actual_coluna=abs(q_actual_linha*q-k);
   printf("q%d k_counter%d\n",q,k_counter);

   for (line_actual=0;line_actual<n_fox;line_actual++){
      for (col_actual=0;col_actual<n_fox;col_actual++){

        colu=k_counter*n_fox+col_actual;

      transition_matrix_temp[POS(line_actual+q_actual_linha,colu)]=transition_matrix_temp[position];
        position++;

        printf(" Process:%d will pass for the final matrix at absolute Position:%d the contents of line %d and column %d of the old matrix\n",k,position,line_actual+q_actual_linha,colu);
     }
  }
  k_counter++;
}
memcpy(transition_matrix,transition_matrix_temp,sizeof(double)*size_graph*size_graph);

}




static int
read_file(const char *fname,double damping,int size,int rank)
{

   FILE *file;
   if (file = fopen(fname, "r"))
   {
      int i,j,total;
//gets size
      char line[1024];
      fgets(line,sizeof(line),file);
      size_graph=atoi(line);
     printf("debug\n");
//      printf("q%d size_graph:%d mod:%d\n",q,size_graph,size_graph%q);

      if (method_flag==1){
         assert((size_graph %q) == 0.0);
         n_fox=size_graph/q;
      }

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

      if(method_flag==1 && rank==0)
         reorder_matrix(size);

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


void zero_matrix(double* submatrix,int size){
   int i,j;
   for(i=0;i<size;i++){
    for(j=0;j<size;j++){
     submatrix[POS(i,j)]=0.0;
  }
}
}

void normal_method(int size,int iterations, int rank){


// number of rows that each process is getting
   const int rows = size_graph / size;

// transmit pagerank vector

   double *submatrix = (double *)malloc(sizeof(double)*size_graph * rows);
   double *partial_pagerank_vector = (double *)malloc(sizeof(double) * rows);


// cycle for iterators


   while(converge==0 && temp_iterator< iterations)
   {


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
MPI_Barrier(MPI_COMM_WORLD);

for(i = 0; i < rows; ++i) {

   double acc = 0.0;
   converge_temp=1.0;
//printf("converge_temp:%d\n",converge_temp);
   for(j = 0; j < size_graph; ++j) {
      acc += submatrix[POS(i, j)] * pagerank_vector[j];
   }

//printf("old:%f new:%f \n ",pagerank_vector[numbers+rows*rank],acc);
   partial_pagerank_vector[i] = acc;
//printf("converge_temp:%d logical %f\n",converge_temp,(float)(fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009));
   converge_temp= converge_temp * (int) (fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009);
//covergence test - working by reduce


// debug messages
//printf("acc:%f item:%d rank:%d before:%f difference:%f boolean :%d converge_temp: %d \n",acc,numbers+rows*rank,rank,pagerank_vector[numbers+rows*rank],fabs(pagerank_vector[numbers+rows*rank]-acc),(fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009),converge_temp);
   numbers++;



   MPI_Bcast(pagerank_vector, size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}
MPI_Barrier(MPI_COMM_WORLD);

/*
MPI_Gather is the inverse of MPI_Scatter.
Instead of spreading elements from one process to many processes, MPI_Gather takes elements
from many processes and gathers them to one single process. This routine is highly useful to
many parallel algorithms, such as parallel sorting and searching.
*/
//used for converge


MPI_Gather(partial_pagerank_vector, rows, MPI_DOUBLE, pagerank_vector, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

MPI_Reduce(&converge_temp, &converge, 1,MPI_INT, MPI_PROD,0,MPI_COMM_WORLD);
if(rank == 0){
   print_pagerank_vector();
   temp_iterator++;
   printf("Iteration nº%d\n",temp_iterator);
   if (converge==1) {
      printf("Values Converged at iteration %d!\n",temp_iterator);
//temp_iterator=iterations;
      converge_temp=1;


   }

}

MPI_Bcast(&temp_iterator, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&converge, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&converge_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);





}

//printf("here");
MPI_Finalize();

// free dynamic memory
free(submatrix);
free(partial_pagerank_vector);
free(pagerank_vector);
if(rank == 0)
   {free(transition_matrix);
   printf("NORMAL METHOD\n");
   }

}

void matrix_multiply(double *A,double *B, double *C){
   int i,j,k;
   for(i=0;i<n_fox;i++){
    for(j=0;j<n_fox;j++){
     //C[i][j]=0;
     for(k=0;k<n_fox;k++){
      C[POS(i,j)]+=A[POS(i,k)]*B[POS(k,j)];
   }
}
}
}

void fox_method(int size,int iterations, int rank, int q){
//MPI_Scatterv
   int temp_iterator=0;
   double acc = 0.0;
   converge_temp=1.0;

   int rankcol, rankline, to, from,root_b;
   int tag=100;




   MPI_Comm_rank(col, &rankcol);
   MPI_Comm_rank(row, &rankline);
   printf("rank:%d col:%d line:%d",rank,rankcol,rankline);

   int* displacements=(int*)malloc(sizeof(int)*q);
   int* scounts=(int*)malloc(sizeof(int)*q);
   int stage;

   const int rows = size_graph / size;

// transmit pagerank vector

   double *submatrixA = (double *)malloc(sizeof(double)*n_fox * n_fox);
   double *submatrixB = (double *)malloc(sizeof(double)*n_fox * n_fox);
   double *submatrixC = (double *)malloc(sizeof(double)*n_fox * n_fox);
   double *submatrixT = (double *)malloc(sizeof(double)*n_fox * n_fox);
   double *submatrix = (double *)malloc(sizeof(double)*n_fox * n_fox);
   int k,l;

   MPI_Scatter(transition_matrix, n_fox*n_fox, MPI_DOUBLE, submatrixA, n_fox*n_fox,MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Scatter(transition_matrix, n_fox*n_fox, MPI_DOUBLE, submatrixB, n_fox*n_fox,MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
   if (rank==0)
   {
      zero_matrix(submatrixC,n_fox);
      zero_matrix(submatrixT,n_fox);
   }

   MPI_Bcast(submatrixC,n_fox*n_fox,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(submatrixT,n_fox*n_fox,MPI_DOUBLE,0,MPI_COMM_WORLD);

   int prev,next;
   while(temp_iterator< iterations)
   {
         
      for (stage=0;stage<q;stage++){
         root_b=(rankline+stage)%q;
         
         prev=(rankline+q-1)%q;
         next=(rankline+1)%q;

            memcpy(submatrixA,submatrixT,sizeof(double)*n_fox*n_fox);
            MPI_Bcast(submatrixT,n_fox*n_fox,MPI_DOUBLE,root_b,row);
            matrix_multiply(submatrixT,submatrixB,submatrixC);

         printf("\n|||||||||||||||||||prev %d next %d actual %d\n", prev,next,rankline);
         MPI_Barrier(row);
         MPI_Sendrecv_replace(submatrixB, n_fox*n_fox, MPI_DOUBLE, prev, tag, next, tag, col, MPI_STATUS_IGNORE);

      }

      printf("submatrixB[0]:%f\n",submatrixB[0]);

      if (rank==0)
      {
         temp_iterator++;
      }

         MPI_Bcast(&temp_iterator, 1, MPI_INT, 0, MPI_COMM_WORLD);
   }

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Gather(submatrixC, n_fox*n_fox, MPI_DOUBLE, transition_matrix, n_fox*n_fox, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (rank==0){
      reorder_matrix_correctly(size_graph);
   }

   MPI_Bcast(transition_matrix, size_graph*size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);


   double *partial_pagerank_vector = (double *)malloc(sizeof(double) * rows);
   MPI_Bcast(pagerank_vector, size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   MPI_Scatter(transition_matrix, rows*size_graph, MPI_DOUBLE, submatrix, rows*size_graph,
   MPI_DOUBLE, 0, MPI_COMM_WORLD);
int numbers=0;


// multiply submatrix
int i, j;
MPI_Barrier(MPI_COMM_WORLD);

for(i = 0; i < rows; ++i) {

   double acc = 0.0;
   converge_temp=1.0;
//printf("converge_temp:%d\n",converge_temp);
   for(j = 0; j < size_graph; ++j) {
      acc += submatrix[POS(i, j)] * pagerank_vector[j];
   }

//printf("old:%f new:%f \n ",pagerank_vector[numbers+rows*rank],acc);
   partial_pagerank_vector[i] = acc;
//printf("converge_temp:%d logical %f\n",converge_temp,(float)(fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009));
   converge_temp= converge_temp * (int) (fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009);
//covergence test - working by reduce


// debug messages
//printf("acc:%f item:%d rank:%d before:%f difference:%f boolean :%d converge_temp: %d \n",acc,numbers+rows*rank,rank,pagerank_vector[numbers+rows*rank],fabs(pagerank_vector[numbers+rows*rank]-acc),(fabs(pagerank_vector[numbers+rows*rank]-acc)<0.0000009),converge_temp);
   numbers++;



   MPI_Bcast(pagerank_vector, size_graph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}
MPI_Barrier(MPI_COMM_WORLD);

/*
MPI_Gather is the inverse of MPI_Scatter.
Instead of spreading elements from one process to many processes, MPI_Gather takes elements
from many processes and gathers them to one single process. This routine is highly useful to
many parallel algorithms, such as parallel sorting and searching.
*/
//used for converge


MPI_Gather(partial_pagerank_vector, rows, MPI_DOUBLE, pagerank_vector, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

MPI_Reduce(&converge_temp, &converge, 1,MPI_INT, MPI_PROD,0,MPI_COMM_WORLD);
if(rank == 0){
   print_pagerank_vector();

}


printf("here");
MPI_Finalize();

// free dynamic memory
free(submatrix);
free(submatrixA);
free(submatrixB);
free(submatrixC);
free(partial_pagerank_vector);
free(pagerank_vector);


if(rank == 0)
   {free(transition_matrix);}
 //MPI_Type_free(submatrix_type);
printf("FOX METHOD n_fox:%d q:%d\n",n_fox,q);
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

   if (strcmp(method,"fox")==0)
      method_flag=1;


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


      if (method_flag==1){

       q= (int)sqrt((double)size);
    }
   
    assert(read_file(file,damping,size,rank)!=0);


//use_dummy_graph(damping);
    print_transition_matrix();
    print_pagerank_vector();

// we ensure that the transition matrix can be evenly divided

    if(method_flag!=1)
       assert(size_graph % size == 0);

    converge=0;
    converge_temp=0;
    temp_iterator=0;
 }

 MPI_Bcast(&size_graph, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&converge, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&converge_temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&temp_iterator, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&n_fox, 1, MPI_INT, 0, MPI_COMM_WORLD);

// non zero processes have not allocated their pagerank vector yet
 if(rank != 0)
   allocate_pagerank_vector();

if (method_flag!=1){

   normal_method(size,iterations,rank);
}

else{
   printf("FOX\n");

   int ndims, reorder, ierr;
   int dim_size[2], periods[2];
   ndims = 2;
   dim_size[0] = q;
   dim_size[1] = q;
   periods[0] = 1;
   periods[1] = 1;
   reorder = 1;

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size,periods, reorder, &new_comm);
   int direction[2] = { 1, 0 };
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Cart_sub(new_comm, direction, &row);
   MPI_Barrier(MPI_COMM_WORLD);
   direction[0] = 0;
   direction[1] = 1;
   MPI_Cart_sub(new_comm, direction, &col);

   MPI_Type_vector(n_fox*n_fox,n_fox,size_graph,MPI_DOUBLE,&submatrix_type);
   MPI_Type_commit (&submatrix_type);
   fox_method(size,iterations,rank,q);

}


return EXIT_SUCCESS;
}

