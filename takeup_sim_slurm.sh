#!/bin/bash

#SBATCH --job-name=takeup_sim    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=48       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=2-00:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=karimn2.0@gmail.com
#SBATCH --output=temp/log/takeup_sim-%j.log
#SBATCH --error=temp/log/takeup_sim-%j.log

module purge
module load rh/devtoolset/8 gdal

Rscript run_stan_dist_sim.R --cmdstanr --include-paths=~/Code/takeup/stan_models --output-path=/tigress/kn6838/takeup --num-sim=12 --num-cores=48 --no-progress-bar --sim-iter=400

