#!/bin/sh
#SBATCH --job-name=ReAct

#### APPROPRIATE FOR A LIKELIHOOD COMPUTATION
#SBATCH --ntasks=1
#SBATCH --time=0-00:10:00

#### APPROPRIATE FOR AN MCMC RUN
###SBATCH --ntasks=32
###SBATCH --time=0-100:00:00

#SBATCH --clusters=baobab
#SBATCH --output=./output/slurm-%J.out
#SBATCH --error=./error/slurm-%J.err

cd /home/bose/cosmosis
bash manual-install-setup
cd /home/bose/cosmosis

### SEQUENTIAL
cosmosis ../ReAct/reactions/tests/test_pyreact/pipeline_lsst.ini
### PARALLEL WITH MPI
##mpirun -n 32 --mca mtl ^psm --mca btl ^openib --oversubscribe cosmosis --mpi  ../ReAct/reactions/tests/test_pyreact/pipeline_lsst.ini
