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
  nparts = np;

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

  elmdist = (idx_t *)malloc(sizeof(idx_t)*(np+1));
  idx_t nepartition = ne/np;
  idx_t neleftover  = ne%np;
  /* if (rank == 0) printf("ne partition = %ld\n", nepartition); */
  /* if (rank == 0) printf("ne left over = %ld\n", neleftover); */
  for (i = 0; i < np; i++) {
    elmdist[i] = nepartition*i;
  }
  elmdist[np] = nepartition*np + neleftover; // fill the last value of elmdist
  int neloc = nepartition;
  if (rank == np - 1) neloc += neleftover;

  //printf("My rank = %d, neloc = %d\n", rank, neloc);

  if (rank == 0) {
  printf("\n");
    for (i=0; i < np + 1; i++){
      printf("%ld ", elmdist[i]);
    }
  printf("\n");
  }


  /* // Allocate eind and eptr arrays */
  eind = (idx_t *)malloc(sizeof(idx_t)*(neloc*4));
  eptr = (idx_t *)malloc(sizeof(idx_t)*neloc +1);

  // Loop over values in file and fill arrays
  int k = 0;
  for (i = 0; i < ne; i++) {
    if (i >= elmdist[rank] && i < elmdist[rank+1]) {
      for (j = 0; j < 4; j++) {
	fscanf(connFile, "%ld", &eind[k*4 + j]);
      }
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
    if (rank == j) printf("eloc array:\n");
    if (rank == j) for (i=0; i < neloc*4; i++) printf("%ld ", eind[i]);
    if (rank == j) printf("\n");
  }
  /* fclose(connFile); */
  /* /\* ----------------------------- END READING FILES ----------------------------------------*\/ */

  /* // Allocate epart and npart arrays based on ne, nn */
  /* epart = (idx_t *)malloc(sizeof(idx_t)*ne); */
  /* npart = (idx_t *)malloc(sizeof(idx_t)*nn); */

  /* int ParMETIS V3 PartMeshKway ( idx_t *elmdist, idx_t *eptr, idx_t *eind, */
  /* 				 idx_t *elmwgt, idx_t *wgtflag, idx_t *numflag, */
  /* 				 idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, */
  /* 				 real_t *tpwgts, real_t *ubvec, */
  /* 				 idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm ) */
    
  /* // Call METIS routine */
  /* METIS_PartMeshDual(&ne, &nn, eptr, eind, */
  /* 		     NULL, NULL, &ncommon, &nparts, */
  /* 		     NULL, NULL, &edgecut, epart, npart ); */

  /* printf("Computed edgecut = %d\n", edgecut); */

  /* free(eptr);  free(eind);  free(epart);  free(npart); */

  //free(elmdist);
  MPI_Finalize();

  return 0;
}
