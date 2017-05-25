#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "parmetis.h"

int main(int argc, char* argv[]) {
  idx_t  nn, ne, edgecut, nparts = 1, count = 0;
  idx_t  *eind, *eptr;
  char   *connFileName, *nodeFileName, *line = NULL;
  FILE   *connFile,  *nodeFile;
  size_t len = 0;
  int    i, j, testing = 0;

  idx_t *elmdist;
  
  MPI_Init(NULL, NULL);
  int np, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  nparts = (idx_t)np;

  // Get command line arguments
  if (argc < 4) {
    if (rank == 0) printf("Usage:\n\t./<exe-name> n ../meshes/<connfile>.conn ../meshes/<nodefile>.xyz \n");
    exit(1);
  } else {
    connFileName = argv[2];
    nodeFileName = argv[3];
  }
  if (argc == 5) testing = atoi(argv[4]);

  // Read in number of nodes
  nodeFile = fopen(nodeFileName, "r");
  getline(&line, &len, nodeFile);
  nn = atol(line);
  nn *= 1;
  fclose(nodeFile);
  
  // Read in number of elements
  connFile = fopen(connFileName, "r");
  fscanf(connFile, "%d", &ne);

  elmdist = (idx_t *)malloc(sizeof(idx_t)*(nparts+1));
  idx_t nepartition = ne/nparts;
  idx_t neleftover  = ne%nparts;

  for (i = 0; i < np; i++) {
    elmdist[i] = nepartition*i;
  }
  elmdist[nparts] = nepartition*nparts + neleftover; // fill the last value of elmdist

  idx_t neloc = nepartition;
  if (rank == nparts - 1) neloc += neleftover; // Give the last node the remainder of elements

  //printf("My rank = %d, neloc = %d\n", rank, neloc);

  if (testing == 1){
    if (rank == 0) {
      printf("elmdist array:\n");
      for (i=0; i < nparts + 1; i++){
	printf("%d ", elmdist[i]);
      }
      printf("\n");
    }
  }

  // Allocate eind and eptr arrays
  eind = (idx_t *)malloc(sizeof(idx_t)*(neloc*4));
  eptr = (idx_t *)malloc(sizeof(idx_t)*neloc +1);

  // Loop over values in file and fill arrays
  int k = 0;
  idx_t tmpval;
  for (i = 0; i < ne; i++) {
    //printf("\nrank %d, elmdist[rank]= %d, i=%d\n", rank, elmdist[rank], i);
    for (j = 0; j < 4; j++) {
      fscanf(connFile, "%d", &tmpval);
      if (i >= (int) elmdist[rank] && i < (int)elmdist[rank+1]) {
	eind[k*4 + j] = tmpval;
	//if (rank == 1) printf("%d ", eind[k*4+j]);
      }
    }
    if (i >= (int) elmdist[rank] && i < (int)elmdist[rank+1]) {
      eptr[k] = count;
      
      count += 4;
      k += 1;
    }
  }

  eptr[neloc] = count; // Fill last value of eind */

  if (testing == 1) {
    for (j = 0; j < np; j++) {
      if(rank == j) printf("\nRank %d\n",j);
      if(rank == j)   printf("N elements = %d\n", neloc);
      if(rank == j)   printf("eptr array:\n");
      if (rank == j) for (i=0;i<neloc+1;i++) printf("%d ",eptr[i]);
      if (rank == j) printf("\n");
      if (rank == j) printf("eind array:\n");
      if (rank == j) for (i=0; i < neloc*4; i++) printf("%d ", eind[i]);
      if (rank == j) printf("\n");
    }
  }

  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon    = 1;
  idx_t ncommonnodes = 3; // Could be 2 ?
  real_t *tpwgts  = (real_t *)malloc(sizeof(real_t)*ncon*nparts);
  real_t *ubvec   = (real_t *)malloc(sizeof(real_t)*ncon);
  idx_t  *options = (idx_t  *)malloc(sizeof(idx_t));
  idx_t  *part    = (idx_t  *)malloc(sizeof(idx_t) *4*neloc);

  
  options[0] = 0;
  ubvec[0] = (real_t) 1.05;
  for (i=0; i < ncon*nparts; i++) tpwgts[i] = 1.0/nparts;

  /* if (rank == 0) for (i=0; i<ncon*nparts; i++) printf("%f ", tpwgts[i]); */
  /* if (rank == 0) printf("\n");   */

  clock_t begin, end;
  double time_spent;
  if( rank == 0) begin = clock();
  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL,
			   &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts,
			   tpwgts, ubvec, options, &edgecut, part, &mpi_comm);
  
  if (rank == 0) end = clock();
  if (rank == 0) time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  
  if (rank == 0) printf("%s time = %f\n", &argv[2][10],time_spent);
  if (rank == 0) printf("%s edgecut = %d\n", &argv[2][10], edgecut);

  
  if (testing == 1) {
    printf("\n");
    if (rank == 0) for (i=0;i<4*neloc;i++) printf("%d ", part[i]);
    if (rank == 1) for (i=0;i<4*neloc;i++) printf("%d ", part[i]);
    printf("\n");
  }

  free(elmdist);  free(eptr);  free(eind);  free(tpwgts); free(ubvec); free(options); free(part);

  MPI_Finalize();

  return 0;
}
