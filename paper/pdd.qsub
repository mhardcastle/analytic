#!/bin/sh -f
#PBS -N analytic
#PBS -m a
#PBS -l nodes=1:ppn=1
#PBS -q main
#PBS -l walltime=12:00:00
#PBS -W group_list=lofar
#PBS -k oe

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: array ID is $PBS_ARRAYID
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo --------------------------------------------

cd /home/mjh/pdd
export PYTHONPATH=/home/mjh/git/analytic
python /home/mjh/git/analytic/paper/pdd-qsub.py $PBS_ARRAYID
