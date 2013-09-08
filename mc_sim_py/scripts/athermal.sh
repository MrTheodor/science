#!/bin/bash

#$ -cwd
#$ -m bes
#$ -M jkrajniak@gmail.com
#$ -V
#$ -pe mpich2 4

MY_HOST="`hostname`"
MY_DATE="`date`"

echo "Running on $MY_HOST at $MY_DATE"
echo "Running environment:" 
env
echo "======================================="
echo "PE_HOSTFILE file:"
cat $PE_HOSTFILE
echo "======================================="
echo JOB_NAME="$JOB_NAME"
echo JOB_ID="$JOB_ID"
echo QUEUE="$QUEUE"
echo SGE_CWD_PATH=$SGE_CWD_PATH
echo PATH=$PATH
echo NSLOTS=$NSLOTS
echo SGE_STDIN_PATH=$SGE_STDIN_PATH
echo SGE_STDOUT_PATH=$SGE_STDOUT_PATH
echo SGE_STDERR_PATH=$SGE_STDERR_PATH
echo SGE_O_HOST=$SGE_O_HOST
echo SGE_O_PATH=$SGE_O_PATH
echo MPICH_HOME=$MPICH_HOME
echo MY_TEST=$MY_TEST
echo "================================================================"

mpiexec -n $NSLOTS python code/intraglobular/src/run.py Experiment1:athermal run
