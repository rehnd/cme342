#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "parmetis.h"

int main(idx_t argc, char* argv[]) {
  idx_t  nn, ne, edgecut, nparts = 1, count = 0, ncommon = 1;
  idx_t  *eind, *eptr, *epart, *npart;
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

  /* ----------------------------------READ FILES FROM COMMANDLINE -----------------------------*/
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
  if (rank == 0) printf("Found # Nodes    = %ld\n", nn);
  fclose(nodeFile);
  
  // Read in number of elements
  connFile = fopen(connFileName, "r");
  fscanf(connFile, "%ld", &ne);
  if (rank == 0) printf("Found # elements = %ld\n", ne);

  elmdist = (idx_t *)malloc(sizeof(idx_t)*(nparts+1));
  idx_t nepartition = ne/nparts;
  idx_t neleftover  = ne%nparts;
  /* if (rank == 0) printf("ne partition = %ld\n", nepartition); */
  /* if (rank == 0) printf("ne left over = %ld\n", neleftover); */
  for (i = 0; i < np; i++) {
    elmdist[i] = nepartition*i;
  }
  elmdist[nparts] = nepartition*nparts + neleftover; // fill the last value of elmdist
  int neloc = nepartition;
  if (rank == nparts - 1) neloc += neleftover; // Give the last node the remainder of elements

  //printf("My rank = %d, neloc = %d\n", rank, neloc);

  /* if (rank == 0) { */
  /*   for (i=0; i < nparts + 1; i++){ */
  /*     printf("%ld ", elmdist[i]); */
  /*   } */
  /*   printf("\n"); */
  /* } */


  /* // Allocate eind and eptr arrays */
  eind = (idx_t *)malloc(sizeof(idx_t)*(neloc*4));
  eptr = (idx_t *)malloc(sizeof(idx_t)*neloc +1);

  idx_t tmpval;
  // Loop over values in file and fill arrays
  int k = 0;
  for (i = 0; i < ne; i++) {
    //printf("\nrank %d, elmdist[rank]= %ld, i=%d\n", rank, elmdist[rank], i);
      for (j = 0; j < 4; j++) {
	fscanf(connFile, "%ld", &tmpval); //&eind[k*4 + j]);
	if (i >= (int) elmdist[rank] && i < (int)elmdist[rank+1]) {
	  eind[k*4 + j] = tmpval;
	  //if (rank == 1) printf("%ld ", eind[k*4+j]);
	}
      }
      if (i >= (int) elmdist[rank] && i < (int)elmdist[rank+1]) {
	eptr[k] = count;
	
	count += 4;
	k += 1;
      }
  }

  eptr[neloc] = count; // Fill last value of eind */

  for (j = 0; j < np; j++) {
    if(rank == j) printf("\nRank %d\n",j);
    if(rank == j)   printf("N elements = %d\n", neloc);
    if(rank == j)   printf("eptr array:\n");
    if (rank == j) for (i=0;i<neloc+1;i++) printf("%ld ",eptr[i]);
    if (rank == j) printf("\n");
    if (rank == j) printf("eind array:\n");
    if (rank == j) for (i=0; i < neloc*4; i++) printf("%ld ", eind[i]);
    if (rank == j) printf("\n");
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
  for (i=0; i<(int)ncon*nparts; i++) tpwgts[i] = (real_t)1.0/(real_t)nparts;
  
  int ret = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL,
				     &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts,
				     tpwgts, ubvec, options, &edgecut, part, &mpi_comm);

  
  if(rank == 0) printf("Computed edgecut = %ld\n", edgecut);

  printf("\n\n");
  if (rank == 0) for (i=0;i<4*neloc;i++) printf("%ld ", part[i]);
  //if (rank == 1) for (i=0;i<4*neloc;i++) printf("%ld ", part[i]);
  printf("\n");
  free(eptr);  free(eind);  free(elmdist);

  MPI_Finalize();

  return 0;
}
