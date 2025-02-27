#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:20:00
#SBATCH -o results/out_ff.log
#SBATCH -e results/out_ff.err

for i in {8,16,32,64,128,256,1024,2048,4096,6144,8192};do 
	for j in {2,4,6,8,10,12,14,16,20,24,28};do 
	./ff $i $j
	done 
done

