#!/bin/bash
#SBATCH -J 4Stripes_Re2.5E4_Outercylinder
#SBATCH -o TC_job.o%j
#SBATCH -t 72:00:00
#SBATCH -N 1 -n 20

module load fftw/3.3.4 intellib intellib/ics-2016 intelmpi/5.1.1 hdf5/1.8.17 fftw/3.3.4

mpirun  -np 20 ./boutnp

