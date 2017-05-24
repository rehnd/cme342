#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "metis.h"

int main(int argc, char* argv[]) {
  idx_t  nn, ne, edgecut, nparts = 1, count = 0, ncommon = 1;
  idx_t  *eind, *eptr, *epart, *npart;
  char   *connFileName, *nodeFileName, *line = NULL;
  FILE   *connFile,  *nodeFile;
  size_t len = 0;
  int    i, j, testing = 0;

  /* ----------------------------------READ FILES FROM COMMANDLINE -----------------------------*/
  // Get command line arguments
  if (argc < 4) {
    printf("Usage:\n\t./<exe-name> n ../meshes/<connfile>.conn ../meshes/<nodefile>.xyz \n");
    exit(1);
  } else {
    nparts = atoi(argv[1]);
    connFileName = argv[2];
    nodeFileName = argv[3];
  }
  if (argc == 5) testing = atoi(argv[4]);

  // Read in number of nodes
  nodeFile = fopen(nodeFileName, "r");
  getline(&line, &len, nodeFile);
  nn = atol(line);
  printf("Found # Nodes    = %d\n", nn);
  fclose(nodeFile);
  
  // Read in number of elements
  connFile = fopen(connFileName, "r");
  fscanf(connFile, "%d", &ne);
  printf("Found # elements = %d\n", ne);

  // Allocate eind and eptr arrays
  eind = (idx_t *)malloc(sizeof(idx_t)*(ne*4));
  eptr = (idx_t *)malloc(sizeof(idx_t)*ne +1);

  // Loop over values in file and fill arrays
  for (i = 0; i < ne; i++) {
    for (j = 0; j < 4; j++) {
      fscanf(connFile, "%d", &eind[i*4 + j]);
      if (testing) printf("%d ", eind[i*4+j]);
    }
    eptr[i] = count;
    if (testing) printf("    %d\n", eind[i]);
    count += 4;
  }
  eptr[i] = count; // Fill last value of eind

  fclose(connFile);
  /* ----------------------------- END READING FILES ----------------------------------------*/

  // Allocate epart and npart arrays based on ne, nn
  epart = (idx_t *)malloc(sizeof(idx_t)*ne);
  npart = (idx_t *)malloc(sizeof(idx_t)*nn);

  // Call METIS routine
  METIS_PartMeshDual(&ne, &nn, eptr, eind,
		     NULL, NULL, &ncommon, &nparts,
		     NULL, NULL, &edgecut, epart, npart );

  printf("Computed edgecut = %d\n", edgecut);

  free(eptr);  free(eind);  free(epart);  free(npart);

  return 0;

}
