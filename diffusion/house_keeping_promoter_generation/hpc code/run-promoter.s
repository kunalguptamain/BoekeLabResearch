#!/bin/bash -e

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=7:00:00
#SBATCH --mem=100GB
#SBATCH --gres=gpu:a100:2
#SBATCH --job-name=promoter_diffusion
#SBATCH --mail-type=END
#SBATCH --mail-user=kg3163@nyu.edu
#SBATCH --output=promoter_diffusion_%A_%a.out
#SBATCH --error=promoter_diffusion_%A_%a.err
#SBATCH --array=0

export PYTHONUNBUFFERED=TRUE

# Define and create a unique scratch directory for this job
mkdir -p "${SLURM_JOB_NAME}_${SLURM_JOBID}" && cd "${SLURM_JOB_NAME}_${SLURM_JOBID}" || exit -1

module purge

# Set TF_GPU_ALLOCATOR environment variable
##export TF_GPU_ALLOCATOR=cuda_malloc_async

# Nvidia GPUs gpu:a100:3 (PyTorch 2.1)
singularity exec --nv \
	    --overlay /scratch/kg3163/promoters/env/overlay-15GB-500K.ext3 \
	    /scratch/work/public/singularity/cuda11.8.86-cudnn8.7-devel-ubuntu22.04.2.sif \
	    /bin/bash -c "source /ext3/env.sh; python ${SLURM_SUBMIT_DIR}/diffusion_promoter.py"


# Nvidia GPUs gpu:rtx8000:1 (TensorFlow 2.13)
#singularity exec --nv \
#	    --overlay /scratch/et1799/Smart_Wearable_Devices/Classifiers/tf2.13cuda11.8.ext3:ro \
#	    /scratch/work/public/singularity/cuda11.8.86-cudnn8.7-devel-ubuntu22.04.2.sif \
#	    /bin/bash -c "source /ext3/env.sh; python ${SLURM_SUBMIT_DIR}/ES_3DCapsNet.py"


# Nvidia GPUs gpu:rtx8000:1
#singularity exec --nv \
#	    --overlay /scratch/et1799/Smart_Wearable_Devices/Classifiers/tf2.4cuda11.0.ext3:ro \
#	    /scratch/work/public/singularity/cuda11.0-cudnn8-devel-ubuntu18.04.sif \
#	    /bin/bash -c "source /ext3/env.sh; python ${SLURM_SUBMIT_DIR}/ES_3DCapsNet.py"

# AMD GPUs gpu:mi50:1
#singularity exec --rocm \
#	    --overlay /scratch/et1799/Smart_Wearable_Devices/Classifiers/tf2.10_amd.ext3:ro \
#	    /scratch/work/public/singularity/rocm5.2.3-ubuntu20.04.5.sif \
#	    /bin/bash -c "source /ext3/env.sh; python ${SLURM_SUBMIT_DIR}/ES_3DCapsNet.py"



