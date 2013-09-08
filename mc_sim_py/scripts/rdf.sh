#!/bin/bash

#$ -cwd
#$ -m bes
#$ -M jkrajniak@gmail.com
#$ -V
#$ -j y
#$ -notify
#$ -pe mpich2 1

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

#INITIAL_CONFIGURATION="`pwd`/17_box_999999.dump"
#INITIAL_CONFIGURATION="`pwd`/experiment/exp5_0/const_t_1_14_32_32_32_5/data/box_900000.dump"
#echo INITIAL_CONFIGURATION="$INITIAL_CONFIGURATION"

mpiexec -n $NSLOTS python ../code/intraglobular/src/tools/calculate_rdf_mpi.py exp8c_xyz_new_2

echo "STOP on `date`"
post_process.sh $JOB_ID $MY_LOG_DIR
