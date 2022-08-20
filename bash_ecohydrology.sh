#!/bin/bash

#SBATCH --job-name=scot   		## Name of the job.
#SBATCH -p free               	## partition/queue name
#SBATCH -N 1         			## run on a single node
#SBATCH -n 1         			## request 12 tasks (12 CPUs)
#SBATCH -t 01:20:00  			## 80 min run time limit

module load R/3.6.2
R CMD BATCH --no-save mycode.R