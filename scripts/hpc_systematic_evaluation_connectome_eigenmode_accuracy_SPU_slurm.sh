#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-16:00:00
#SBATCH --mem=60G
#SBATCH --partition=cascade
#SBATCH -o /data/gpfs/projects/punim0695/fMRI/slurm/eigenmodes/accuracy_%j.out

module load OpenSSL/1.1 foss/2022a Python/3.10.4
cd /data/gpfs/projects/punim0695/fMRI
source new_venv/bin/activate

# some colors for fancy logging :D
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# read inputs
density=$1
binary=$2
regression=$3
tractography=$4
fwhm=$5
global_local_epsilon=$6

echo -e "${GREEN}[INFO]${NC} `date`: Starting..."
echo -e "${GREEN}[INFO]${NC} `date`: Parameters: \n\t${density} \n\t${binary} \n\t${regression} \n\t${tractography} \n\t${fwhm} \n\t${global_local_epsilon}"

python3 hpc_systematic_evaluation_connectome_eigenmode_accuracy.py "${density}" "${binary}" "${regression}" "${tractography}" "${fwhm}" "${global_local_epsilon}"

echo -e "${GREEN}[INFO]${NC} `date`: Script finished."
