/* 
 * invertastic.c
 * Copyright 2015 The University of Edinburgh

 * Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at 

 * http://www.apache.org/licenses/LICENSE-2.0 

 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
 * See the License for the specific language governing permissions and 
 * limitations under the License.  
 */


/* SEE THE README FOR MORE INFO */

 #include "invertastic.h" 

 static char* usage=" \ 
  \n \ 
  invertastic --size <size> \n \ 
        [--check] \n \ 
        [--input <full path to input file>] \n \ 
        [--output <full path to output file>] \n \ 
  \n \ 
 \n where <size> is N for an NxN square matrix\n"; 


int main(int argc, char *argv[])
{


  int i, j, k;
  int matSize=0;
  int info;

  // parallelisation
  int nPRow=0,nPCol=0,myPRow=0, myPCol=0;
  int myRank=0, numTask=0, context=0, master=0,myId=0,nProc=0;
  unsigned long int nLocRow=0,nLocCol=0;
  int desc[9];

  //timing
  double time0, time1, time2;

  // I/O
  FILE *fp;
  char* inputFile_ptr;
  char* outputFile_ptr;
  char* inputFile;
  char* outputFile;
  int inputFlag=0;
  int outputFlag=0;

  //correctness checking
  int checkFlag=0;



 // ---General Parallelisation setup----

   // initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&numTask);
    
  // query BLACS for the number of tasks, and the id of this task
    Cblacs_pinfo( &myId, &nProc );

  // set master rank to 1
  if (myRank==0) master=1;

  if(master) printf("=================================================\n");
  if(master) printf("             *** Invertastic ***\n");
  if(master) printf("  (C) The University of Edinburgh (2015)\n");
  if(master) printf("  Please see the README for more information \n");
  if(master) printf("=================================================\n\n");



  // get a context handle from BLACS
  Cblacs_get( 0, 0, &context );

  //initialise BLACS to a 2D regular grid
  nPRow=sqrt(numTask);

  while (numTask%nPRow)nPRow++;

  nPCol=numTask/nPRow;

  if (numTask > 1){

  if (nPRow*nPCol != numTask || nPRow == numTask){
    if (master) printf("Error:could not create a rectangular process grid \n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  }
  

  if(master) printf("Running with %d tasks, decomposed in a %dx%d grid\n\n",numTask,nPRow,nPCol);

  //processor decomposition is row major
  Cblacs_gridinit( &context, "Row-major", nPRow, nPCol );
  //Cblacs_gridinit( &context, "Column-major", nPRow, nPCol );


 // query BLACS for the processor grid information
  Cblacs_pcoord( context, myId, &myPRow, &myPCol);

  // ------End General Parallelisation setup----------



  // ------Parse command line arguments-------

  //master task gets details and broadcasts to other tasks
  if (master){

    for (i=1;i<argc;i++){
 

     if (!strcmp(argv[i],"--size")){
	i++;
	if(i<argc)
	  matSize=atoi(argv[i]);
	printf("Matrix Size is %dx%d \n\n",matSize,matSize);
	continue;
      }

     if (!strcmp(argv[i],"--input")){
	i++;
	if(i<argc)
	  inputFile_ptr=argv[i];
	printf("Input file is %s \n\n",inputFile_ptr);
	inputFlag=1;
	continue;
      }


     if (!strcmp(argv[i],"--output")){
	i++;
	if(i<argc)
	  outputFile_ptr=argv[i];
	printf("Output file is %s \n\n",outputFile_ptr);
	outputFlag=1;
	continue;
      }

     if (!strcmp(argv[i],"--check")){
	printf("--check specified: resulting inverse will be multiplied by original and compared with identity\n\n");
	checkFlag=1;
	continue;
      }

     printf("Error. Argument %s not recognised. Usage:\n %s\n",argv[i],usage);
     MPI_Abort(MPI_COMM_WORLD,0);

      
    }
        
    if (matSize==0){
      printf("Error. Usage:\n %s\n",usage);
      MPI_Abort(MPI_COMM_WORLD,0);
    }


  }

  MPI_Bcast( &matSize, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD ); 

  MPI_Bcast( &inputFlag, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD ); 
  
  MPI_Bcast( &outputFlag, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD ); 

  MPI_Bcast( &checkFlag, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD ); 


  int filenameLength=0;

  if (inputFlag){

    if (master) filenameLength=strlen(inputFile_ptr)+1;
    MPI_Bcast( &filenameLength, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD );

    inputFile=malloc(filenameLength*sizeof(char));
    if (master) strcpy(inputFile,inputFile_ptr);
    MPI_Bcast( inputFile, filenameLength, MPI_CHAR, 0,  MPI_COMM_WORLD );

  }


  if (outputFlag){

    if (master) filenameLength=strlen(outputFile_ptr)+1;
    MPI_Bcast( &filenameLength, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD );

    outputFile=malloc(filenameLength*sizeof(char));
    if (master) strcpy(outputFile,outputFile_ptr);
    MPI_Bcast( outputFile, filenameLength, MPI_CHAR, 0,  MPI_COMM_WORLD );

  }


  MPI_Barrier(MPI_COMM_WORLD);

  // ------End Parse command line arguments-------


  //-----Set up distributed matrix---
  
 // get the number of rows and columns for local matrix
  nLocRow = numroc_( &matSize, &NBRow, &myPRow, &ZERO, &nPRow );
  nLocCol = numroc_( &matSize, &NBCol, &myPCol, &ZERO, &nPCol );


  // set up the matrix descriptor
  descinit_(desc, &matSize,   &matSize,   &NBRow,  &NBCol,
 	    &ZERO, &ZERO, &context, &nLocRow,  &info);
   if(info!=0){
     if (master) printf("Error: desc init - return value %d \n",info);
     MPI_Abort(MPI_COMM_WORLD,1);
   }

  //-----End set up distributed matrix---


   //----set up the MPI-IO------


   //create block cyclic distributed array datatype
   //we use fortran (column major) ordering for data for ScaLAPACK compatibility
   int dims[2] = {matSize,matSize};
   int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
   int dargs[2] = {NBRow, NBCol};
   int pdims[2]={nPRow, nPCol};
   MPI_Datatype darray;

   MPI_Type_create_darray(numTask, myRank, 2, dims, distribs, dargs, 
			  pdims, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,
			  &darray);
   MPI_Type_commit(&darray);
   
   int locsize;
   
   //get size in bytes of type
   MPI_Type_size(darray, &locsize);

   //and corresponding number of elements
   int nelements = locsize/8;
   
   //get lower bound and extent
   MPI_Aint lb, locextent;
   MPI_Type_get_extent(darray, &lb, &locextent);
   
   //declare handles for later usage
   MPI_File fileid;
   MPI_Status status;
   
   MPI_Barrier(MPI_COMM_WORLD);
   
   //----end set up the MPI-IO------
    

   // ----allocate memory for local part of matrix---

   //if(master) printf("local malloc size: %d \n",nLocRow*nLocCol*sizeof(double));
   
   double *matrix = (double*) calloc (nLocRow*nLocCol,sizeof(double));
   double *matrix_save = (double*) malloc (nLocRow*nLocCol*sizeof(double));
   int *ipiv = (int*) calloc( nLocRow, sizeof(int) );
   double *work = (double*) calloc( 10, sizeof(double));
   int *iwork = (int*) calloc( 10, sizeof(int) );
   
   if ( (!matrix) || (!matrix_save) || (!ipiv) || (!work) || (!iwork) ) {
     if (master) printf("Error. Memory allocation failed.\n");
     MPI_Abort(MPI_COMM_WORLD,1);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   
   // ----end allocate memory for local part of matrix---
   


  // ----populate distributed matrix---
  if(!inputFlag){ //no input file has been specified, so create random SPD matrix
    
    if(master) printf("No input file specified. ");
    if(master) printf("Creating Random %dx%d matrix..\n",matSize,matSize);
    
    for (i = 0; i < nLocRow; i++) 
      for (j = 0; j < nLocCol; j++) 
	matrix[i*nLocCol+j] = ((double) rand())/RAND_MAX;            
    
    if(master) printf("...Done\n\n");
    
    //copy to arrayd_save
    memcpy(matrix_save,matrix,nLocCol*nLocRow*sizeof(double));
    
        
    if(master) printf("Making matrix symmetric positive definite\n");
    //  A=A+A'
    pdgeadd(&trans, &matSize, &matSize, &dONE, matrix_save, &ONE, &ONE, desc, &dONE, matrix, &ONE, &ONE,desc);
    
    
    //A=A+N*eye(N)
    for (i = 0; i < nLocCol; i++) {
      //get global col index
      int globCol = indxl2g(i,NBCol,myPCol,nPCol); 
      for (j = 0; j < nLocRow; j++){
	//get global row index
	int globRow = indxl2g(j,NBRow,myPRow,nPRow); 
	
	if (globRow==globCol)
	  matrix[i* nLocRow+j]  += (double) matSize; //matrix has fortran ordering           
	
	
      }
    }

  MPI_Barrier(MPI_COMM_WORLD);    
  if(master) printf("...Done\n\n");
  }

  else { //read the specified input file 

    if(master) printf("Reading file %s...\n",inputFile);
      

    MPI_Barrier(MPI_COMM_WORLD);      
    time0 = omp_get_wtime();
      

    MPI_File_open(MPI_COMM_WORLD, inputFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fileid);

    if (!fileid){
     if (master) printf("Error: cannot open file \n");
     MPI_Abort(MPI_COMM_WORLD,1);
    }


    if(numTask==1){
      
      MPI_File_read_all(fileid,matrix,matSize*matSize,MPI_DOUBLE_PRECISION,&status);
      
    }
    else{
      
      int disp=0;//bytes to skip
      MPI_File_set_view(fileid, disp, MPI_DOUBLE_PRECISION, darray, "native", MPI_INFO_NULL);
      MPI_File_read_all(fileid,matrix,nelements,MPI_DOUBLE_PRECISION,&status);
      
    }

    MPI_File_close(&fileid);
      

    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    time1 = omp_get_wtime();
    if(master) printf("Reading File: %1.5f seconds\n",time1-time0);  
    if(master) printf("...Done\n\n");
   


    MPI_Barrier(MPI_COMM_WORLD);
  

    
  }


  // ----end populate distributed matrix---




  if (checkFlag){
    MPI_Barrier(MPI_COMM_WORLD);
  //take a copy of matrix for later checking 
  memcpy(matrix_save,matrix,nLocRow*nLocCol*sizeof(double));
  }


  // ----Invert Matrix ---

  if(master) printf("\n\nPerforming Matrix Inversion: double precision \n");

  time0 = omp_get_wtime();

  //Replace matrix with it's LU Decomposition
  pdgetrf_( &matSize, &matSize, matrix, &ONE, &ONE, desc, ipiv, &info );
  if(info!=0){
     if (master) printf("Error: Cholesky factorisation - return value %d \n",info);
     MPI_Abort(MPI_COMM_WORLD,1);
  }

  time1 = omp_get_wtime();
  double time_save=time1-time0;
  if(master) printf("pdpotrf: %1.5f seconds\n",time1-time0);

  time0 = omp_get_wtime();
  //Replace (LU decomposition of) matrix with it's inverse
  printf(">>>> %d %lu  %lu n", matSize, nLocRow, nLocCol);
  int negOne = -1; // c is stoopid.
  pdgetri_( &matSize, matrix, &ONE,&ONE, desc, ipiv, work, &negOne, iwork, &negOne, &info );
  int lwork = (int ) work[0];
  int liwork = (int) iwork[0];
  double *workq = (double*) calloc( lwork, sizeof(double));
  int *iworkq = (int*) calloc(liwork , sizeof(int) );
  
   if (  (!workq) || (!iworkq) ) {
      if (master) printf("Error. Memory allocation failed.\n");
      MPI_Abort(MPI_COMM_WORLD,1);
   }
   pdgetri_( &matSize, matrix, &ONE,&ONE, desc, ipiv, workq, &lwork, iworkq, &liwork, &info );



   if(info!=0){
     if (master) printf("Error: Inverse calculation - return value %d \n",info);
     MPI_Abort(MPI_COMM_WORLD,1);
   }

   time1 = omp_get_wtime();
   
   if(master) printf("pdpotri: %1.5f seconds\n",time1-time0);

   if(master) printf("Total inversion time: %1.5f seconds\n\n",
		     time_save+(time1-time0));


  // ----End Invert Matrix ---

  
  
  // ----check for correctness ---
  
  if (checkFlag){
    
    
    double *eyecalc = (double*) malloc (nLocRow*nLocCol*sizeof(double));

    if (!eyecalc) {
      printf("Error. Memory allocation failed.\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    

    MPI_Barrier(MPI_COMM_WORLD);
    // BLAS matrix multiplication
    if(master) printf("\n\nChecking result: multiplying original matrix by its inverse...\n");
    time0 = omp_get_wtime();    
    
    char side='L';
    pdsymm (&side,&uplo,&matSize, &matSize,
	    &dONE, matrix,&ONE,&ONE, desc, matrix_save,
	    &ONE,&ONE,desc,&dZERO, eyecalc ,&ONE,&ONE,desc);
    
    time1 = omp_get_wtime();
    if(master) printf("Matrix multiplication: %1.5f seconds\n",time1-time0);  
    if(master) printf("...Done\n\n");
    
 
    double tmpdiff;
    double diff=0.;
    

    //loop over local matrix
    for (i = 0; i < nLocCol; i++) {
      //get global col index
      int globCol = indxl2g(i,NBCol,myPCol,nPCol); 
      for (j = 0; j < nLocRow; j++){
      //get global row index
	int globRow = indxl2g(j,NBRow,myPRow,nPRow); 
	
	//Check against identity matrix. 

	if (globRow==globCol)
	  //on diagonal - should be 1.
	  tmpdiff = fabs(eyecalc[i*nLocRow+j] - 1.); 
	else
	  //off diagonal - should be 0.
	  tmpdiff = fabs(eyecalc[i*nLocRow+j] - 0.);
	if (tmpdiff > diff) diff = tmpdiff;      
	
	

      }
    }
    
    if(master) printf("Max element-wise difference of resulting matrix from I is %1.3e\n\n",diff);
    
    free(eyecalc); 
    
  }
  
  // ----end check for correctness ---


  // ----Write matrix to file ---
  if(outputFlag){

    MPI_Barrier(MPI_COMM_WORLD);


    // At the moment, only the lower triangle is populated: populate full matrix

    // Get transpose of matrix
    pdgeadd(&trans, &matSize, &matSize, &dONE, matrix, &ONE, &ONE, desc, &dZERO, matrix_save, &ONE, &ONE,desc);

    // populate upper triangle from transpose 
    for (i = 0; i < nLocCol; i++) {
      //get global col index
      int globCol = indxl2g(i,NBCol,myPCol,nPCol); 
      for (j = 0; j < nLocRow; j++){
	//get global row index
	int globRow = indxl2g(j,NBRow,myPRow,nPRow); 

	if (globRow<globCol)
	  matrix[i* nLocRow+j]=matrix_save[i* nLocRow+j];

      }
    }


    MPI_Barrier(MPI_COMM_WORLD);

    if(master) printf("Writing file %s\n",outputFile);
    time0 = omp_get_wtime();
            

    
    if (master) {
      fp = fopen ( outputFile, "wb" ) ;
      fclose(fp);
      if (!fp){
	printf("Error: cannot open file \n");
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
      
    MPI_File_open(MPI_COMM_WORLD, outputFile, MPI_MODE_WRONLY, MPI_INFO_NULL, &fileid);


    if(numTask==1){

      MPI_File_write_all(fileid,matrix,matSize*matSize,MPI_DOUBLE_PRECISION,&status);

    }
    else{
      
     
     int disp=0;//bytes to skip
     MPI_File_set_view(fileid, disp, MPI_DOUBLE_PRECISION, darray, "native", MPI_INFO_NULL);
     MPI_File_write_all(fileid,matrix,nelements,MPI_DOUBLE_PRECISION,&status);
      
    }
 
    MPI_File_close(&fileid);
     
   
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    time1 = omp_get_wtime();
    if(master) printf("Writing File: %1.5f seconds\n",time1-time0);  
    if(master) printf("...Done\n\n");
    
  }
  // ----End write matrix to file ---


  // ----tidy up---

  free(matrix);
  free(matrix_save);
  
  if(master) printf("Invertastic completed.\n"); 
  if(master) printf("======================\n");


  
  MPI_Finalize();

  // ----end tidy up---
  return 0;
}


