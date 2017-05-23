#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "metis.h"

//void readfiles(idx_t argc, char* argv[], idx_t &nn, idx_t &ne, idx_t &nparts, idx_t *eptr, idx_t *eind);

int main(idx_t argc, char* argv[]) {
  idx_t nparts = 1;
  idx_t nn;
  idx_t ne;
  idx_t objval;

  bool testing = 0;
  
  /* ----------------------------------------------------------------------------------------*/
  // Read in files from commandline
  //readfiles(argc, argv, nn, ne, nparts, eptr, eind);
  char *connFileName;
  char *nodeFileName;
  FILE *connFile;
  FILE *nodeFile;
  char *line = NULL;
  size_t len = 0;

  if (argc < 4) {
    printf("Usage:\n\t./<exe-name> n ../meshes/<connfile>.conn ../meshes/<nodefile>.xyz \n");
    exit(1);
  } else {
    nparts = atoi(argv[1]);
    connFileName = argv[2];
    nodeFileName = argv[3];
  }
  if (argc == 5) {
    testing = argv[4];
  }

  // Get the number of nodes
  nodeFile = fopen(nodeFileName, "r");
  getline(&line, &len, nodeFile);
  nn = atol(line);
  printf("Found # Nodes = %d\n", nn);
  fclose(nodeFile);
  
  // Get number of elements
  connFile = fopen(connFileName, "r");
  fscanf(connFile, "%d", &ne);
  printf("Found # elements = %d\n", ne);
 
  idx_t *eind = (idx_t *)malloc(sizeof(idx_t)*(ne+1));
  idx_t *eptr = (idx_t *)malloc(sizeof(idx_t)*ne*4);
  /* idx_t eind[ne+1]; */
  /* idx_t eptr[ne*4]; */
  
  int i, j;
  idx_t count = 0;
  for (i = 0; i < ne; i++) {
    for (j = 0; j < 4; j++) {
      fscanf(connFile, "%d", &eptr[i*4 + j]);
      if (testing) printf("%d ", eptr[i*4+j]);
    }
    eind[i] = count;
    if (testing) printf("    %d\n", eind[i]);
    count += 4;
  }
  eind[i] = count; // Fill last value of eind
  

  /* =============== TESTING ================ */
  if (testing) {
    printf("eptr array\n");
    for (i = 0; i < ne*4; i++) {
      printf("%d ", eptr[i ]);
    }
    printf("\n");
    printf("eind array\n");
    for (i = 0; i < ne + 1; i++) {
      printf("%d ", eind[i]);
    }
    printf("\n");
  }
  /* =============== END TESTING ================ */
  
  fclose(connFile);
  /* DONE READING FILES */
  /* ----------------------------------------------------------------------------------------*/

  idx_t *epart = (idx_t *)malloc(sizeof(idx_t)*ne);
  idx_t *npart = (idx_t *)malloc(sizeof(idx_t)*nn);
  /* idx_t epart[ne]; */
  /* idx_t npart[nn]; */

  int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  //METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, &objval, epart, npart );
  idx_t ncommon = 1;
  METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, options, &objval, epart, npart );
  printf("objval=%d", objval);

  free(eptr);
  free(eind);
  free(epart);
  free(npart);

  return 0;

}



/* void readfiles(idx_t argc, char* argv[], idx_t &nn, idx_t &ne, idx_t &nparts, idx_t *eptr, idx_t *eind) { */
/*   /\*  */
/*      Takes the commandline arguments for mesh files and # of */
/*      partitions, reads in the # nodes and fills eptr and eind */
/*      according to the input files. */
/*   *\/ */
/*   char *connFileName; */
/*   char *nodeFileName; */
/*   FILE *connFile; */
/*   FILE *nodeFile; */
/*   char *line = NULL; */
/*   size_t len = 0; */

/*   if (argc != 4) { */
/*     printf("Usage:\n\t./<exe-name> n connFileName nodeFileName \n"); */
/*     exit(1); */
/*   } else { */
/*     nparts = atol(argv[1]); */
/*     connFileName = argv[2]; */
/*     nodeFileName = argv[3]; */
/*   } */

/*   // Get the number of nodes */
/*   nodeFile = fopen(nodeFileName, "r"); */
/*   getline(&line, &len, nodeFile); */
/*   nn = atol(line); */
/*   printf("Found # Nodes = %d\n", nn); */
/*   fclose(nodeFile); */
  
/*   // Get number of elements */
/*   connFile = fopen(connFileName, "r"); */
/*   fscanf(connFile, "%d", &ne); */
/*   printf("Found # elements = %d\n", ne); */
  
/*   eind = (idx_t *)malloc(sizeof(idx_t)*(ne+1)); */
/*   eptr = (idx_t *)malloc(sizeof(idx_t)*ne*4); */

/*   int i, j; */
/*   idx_t count = 0; */
/*   for (i = 0; i < ne; i++) { */
/*     for (j = 0; j < 4; j++) { */
/*       fscanf(connFile, "%d", &eptr[i*4 + j]); */
/*       //printf("%d ", eptr[i*4+j]); */
/*     } */
/*     eind[i] = count; */
/*     //printf("    %d\n", eind[i]); */
/*     count += 4; */
/*   } */
/*   eind[i] = count; // Fill last value of eind */

/*   /\* for (i = ne*4 - 10; i < ne*4; i++) { *\/ */
/*   /\*   printf("%d ", eptr[i ]); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */
/*   for (i = 0; i < 4; i++) { */
/*     printf("%d ", eind[i]); */
/*   } */
/*   printf("\n"); */

/*   fclose(connFile); */

/* } */
