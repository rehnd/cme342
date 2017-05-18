
  //! Open and write the Tecplot file. Note that this
  //! format is only tested for Tecplot 360 2010.
  //! The processor number that owns each element
  //! will be written to the cell-center.
  
  ofstream sol_file;
  sol_file.precision(15);
  sol_file.open(filename.c_str(), ios::out);
  sol_file << "TITLE = \"Visualization of METIS partitioning\"" << endl;
  sol_file << "VARIABLES = \"x\", \"y\", \"z\", \"epart\"" << endl;
  sol_file << "ZONE NODES=" << nn <<", ELEMENTS="<< ne;
  sol_file <<", DATAPACKING=BLOCK" << ", ZONETYPE=FEBRICK";
  sol_file <<", VARLOCATION=([4]=CELLCENTERED)" << endl;
  
  //! Write coordinates and mesh processor number.
  
  for (int iNode = 0; iNode < nn; iNode++ ) {
    sol_file << coords[iNode][0] << " ";
    if ((iNode+1)%100 == 0) sol_file << endl;
  }
  sol_file << endl;
  
  for (int iNode = 0; iNode < nn; iNode++ ) {
    sol_file << coords[iNode][1] << " ";
    if ((iNode+1)%100 == 0) sol_file << endl;
  }
  sol_file << endl;
  
  for (int iNode = 0; iNode < nn; iNode++ ) {
    sol_file << coords[iNode][2] << " ";
    if ((iNode+1)%100 == 0) sol_file << endl;
  }
  sol_file << endl;
  
  for (int iElem = 0; iElem < ne; iElem++ ) {
    sol_file << epart[iElem] << " ";
    if ((iElem+1)%100 == 0) sol_file << endl;
  }
  sol_file << endl;
  
  //! Write the element connectivity.
  
  iptr = 0;
  for (int iElem = 0; iElem < ne; iElem++) {
    iptr = iElem*4;
    sol_file <<
    eind[iptr+0]+1 <<" "<< eind[iptr+1]+1 <<" "<<
    eind[iptr+2]+1 <<" "<< eind[iptr+2]+1 <<" "<<
    eind[iptr+3]+1 <<" "<< eind[iptr+3]+1 <<" "<<
    eind[iptr+3]+1 <<" "<< eind[iptr+3]+1 << endl;
  }
  
  //! Close the file and free node memory.
  
  sol_file.close();
