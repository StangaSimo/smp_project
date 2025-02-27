#!/bin/bash

#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --time=00:30:00
#SBATCH -o results/out_omp.log
#SBATCH -e results/out_omp.err

for i in {8,16,32,64,128,256,1024,2048,4096,6144,8192};do 
	mpirun -n 6 omp $i 
done

