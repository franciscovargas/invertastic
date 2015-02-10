#!/bin/bash --login
# 
# run_archer.sh
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
#PBS -N invertastic
#
# specify the number of nodes
#PBS -l select=1
#
#PBS -l walltime=1:00:00
#
# Replace budget code below with your project code (e.g. t01)
#PBS -A z01-cse

#specify N for NxN square matrix

N=1024


# specify a temporary workdir
WORKDIR=/work/z01/z01/alang/INVERTASTIC_SCRATCH


rm -rf $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR


#automatically get the number of cores (24 per node)
NP=`qstat -f $PBS_JOBID | awk '/resources_used.ncpus/ {print $3}'`

#run the code
aprun -n $NP  $PBS_O_WORKDIR/src/invertastic \
    --size $N \
    --check 


#    --input /work/z01/z01/alang/stripedir/spd${N}inv.dat \
#    --output /work/z01/z01/alang/stripedir/spd${N}invinv.dat 


