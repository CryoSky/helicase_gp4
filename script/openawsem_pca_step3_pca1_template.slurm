#!/bin/bash
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --job-name=check
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=sj52@rice.edu
#SBATCH --mail-type=FAIL
#SBATCH --exclude=bc8u23n8
#SBATCH -o trial.out
#SBATCH -e trial.err


date
echo My job ran on: $SLURM_NODELIST
echo My jobid was $SLURM_JOBID
echo My submission directory was $SLURM_SUBMIT_DIR

echo `which python`

#for i in $(seq 1 50);do tail -n 1041700 output_${i}.pdb >> output_tail_total.pdb;done
#python /home/sj52/script/pca/pca_from_pdb_var_length.py output_tail_total.pdb
#python /home/sj52/script/pca/calculate_motion_from_nmd.py pca.nmd 5
#python /home/sj52/script/pca/calculate_pca_weight_from_nmd.py pca.nmd test.png

for i in $(seq INT `echo "INT+0" | bc`);do python /home/sj52/script/pca/pca_projection_back.py /scratch/sj52/helicase_project/20210813_helicase_50bp_dna_fixed_globalpara/efatp/run_job_1/pca.nmd /scratch/sj52/helicase_project/20210813_helicase_50bp_dna_fixed_globalpara/afatp/run_job_1/output_${i}.pdb 1 pca1_output_${i}.txt>> pca1_output_${i}.log;done
