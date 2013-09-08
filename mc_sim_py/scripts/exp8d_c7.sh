#!/bin/bash

#$ -cwd
#$ -m bes
#$ -M jkrajniak@gmail.com
#$ -V
#$ -j y
#$ -notify
#$ -pe mpich2 4

MY_HOST="`hostname`"
MY_DATE="`date`"

MY_LOG_DIR="${JOB_ID}_log"
mkdir $MY_LOG_DIR

export SGE_STDOUT_PATH=$SGE_CWD_PATH/$MY_LOG_DIR/${JOB_NAME}.o${JOB_ID}
export SGE_STDERR_PATH=$SGE_CWD_PATH/$MY_LOG_DIR/${JOB_NAME}.e${JOB_ID}

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

#INITIAL_CONFIGURATION="`pwd`/box_2000.dump"
#INITIAL_CONFIGURATION="`pwd`/experiment/exp8_0/ath_128_128_128_0/data/box_999999.dump"
#INITIAL_CONFIGURATION="`pwd`/experiment/exp8c_3/m5_240_240_240_1/data/box_4999.dump"
echo INITIAL_CONFIGURATION="$INITIAL_CONFIGURATION"


mpiexec -n $NSLOTS python ../code/intraglobular/src/run.py exp8:Experiment8d.c7 run

echo "STOP on `date`"
post_process.sh $JOB_ID $MY_LOG_DIR
