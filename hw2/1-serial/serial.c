#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "metis.h"

void readfiles(int argc, char* argv[], int &nn, int &ne, int &nparts, int *eptr, int *eind);

int main(int argc, char* argv[]) {
  int nparts;
  int nn;
  int ne;
  int *eptr;
  int *eind;
  int objval;
  int *epart;
  int *npart;

  readfiles(argc, argv, nn, ne, nparts, eptr, eind);

  epart = (int *)malloc(sizeof(int)*ne);
  npart = (int *)malloc(sizeof(int)*nn);

  METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL, &objval, epart, npart );

  return 0;

}



void readfiles(int argc, char* argv[], int &nn, int &ne, int &nparts, int *eptr, int *eind) {
  /* 
     Takes the commandline arguments for mesh files and # of
     partitions, reads in the # nodes and fills eptr and eind
     according to the input files.
  */
  char *connFileName;
  char *nodeFileName;
  FILE *connFile;
  FILE *nodeFile;
  char *line = NULL;
  size_t len = 0;

  if (argc != 4) {
    printf("Usage:\n\t./<exe-name> n connFileName nodeFileName \n");
    exit(1);
  } else {
    nparts = atoi(argv[1]);
    connFileName = argv[2];
    nodeFileName = argv[3];
  }

  // Get the number of nodes
  nodeFile = fopen(nodeFileName, "r");
  getline(&line, &len, nodeFile);
  nn = atoi(line);
  printf("Found # Nodes = %d\n", nn);
  fclose(nodeFile);
  
  // Get number of elements
  connFile = fopen(connFileName, "r");
  fscanf(connFile, "%d", &ne);
  printf("Found # elements = %d\n", ne);
  
  eptr = (int *)malloc(sizeof(int)*ne*4);
  eind = (int *)malloc(sizeof(int) * (ne+1));

  int i, j;
  int count = 0;
  for (i = 0; i < ne; i++) {
    for (j = 0; j < 4; j++) {
      fscanf(connFile, "%d", &eptr[i + j]);
    }
    eind[i] = count;
    count += 4;
  }
  eind[i] = count; // Fill last value of eind

  fclose(connFile);

}
