# 
# Makefile
# Copyright 2015 The University of Edinburgh
#
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at 
#
# http://www.apache.org/licenses/LICENSE-2.0 
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License.  
#

CC=mpiicc
CFLAGS=-O3 -I$(MKLROOT)/include -openmp 


LFLAGS=  $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a  -lm 


MF= Makefile

EXE= invertastic

SRC= invertastic.c 

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o: 
	$(CC) $(CFLAGS) -c $< -o $(OBJ)


all:	$(EXE)

$(EXE):	$(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) invertastic.h 


clean:
	rm -f $(OBJ)  $(EXE) core

