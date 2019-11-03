#!/bin/sh
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=1:jobtype=gpu:ngpus=1
#PBS -l walltime=00:20:00

cd $PBS_O_WORKDIR

readonly HOME_DIR=/home/users/bip
readonly PROJ_DIR=${HOME_DIR}/simulation4
readonly PROGRAMS_DIRNAME=programs

source /home/users/bip/.miniconda3/bin/activate
python ${PROJ_DIR}/${PROGRAMS_DIRNAME}/equil_only_pep.py

