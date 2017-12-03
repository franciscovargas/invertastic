/* 
 * invertastic.h
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

//#include <mpi.h>
#include <stdio.h>
//#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mkl_lapack.h>
#include <mkl_cblas.h>
#include <string.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>
#include <mpi.h>

#include <sys/time.h>


int ZERO=0;
int ONE=1;
double dZERO=0.;
double dONE=1.;
char trans='T';
char notrans='N';
char uplo='L';

// block cyclic distribution block size
int NBRow=128;
int NBCol=128;

// ScaLAPACK local to global index mapping (for either of the 2 dimensions) 
static inline int indxl2g(int indxloc, int nb, int iproc, int nprocs){
  return nprocs*nb*(indxloc / nb) + 
    indxloc%nb + ((nprocs + iproc - 0) % nprocs)*nb;
}

// ScaLAPACK global index to processor mapping (for either of the 2 dimensions) 
static inline int indxg2p( int indxglob, int nb,  int nprocs ){
  return (indxglob / nb)% nprocs ;
}

