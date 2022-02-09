#!/bin/bash
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --job-name=heli_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=sj52@rice.edu
#SBATCH --mail-type=FAIL
#SBATCH -o trial.out
#SBATCH -e trial.err


date
echo My job ran on: $SLURM_NODELIST
echo My jobid was $SLURM_JOBID
echo My submission directory was $SLURM_SUBMIT_DIR

echo `which python`
module purge
module load GCC/9.3.0 CUDA/11.0.182


for i in $(seq 1 10);do cd run_job_1; j=$(printf "%03d" "$i"); cp ../full_initial_structure/model_${j}.pdb clean.pdb;k=$(echo "1-0.02*${i}" | bc | awk '{printf "%.2f\n", $0}'); python protein_dna_ligand_run_qdiff.py ${k};mv output.pdb output_${i}.pdb;mv output.dcd output_${i}.dcd;cd ..;done
#for i in $(seq 27 30);do cd run_job_1; j=$(printf "%03d" "$i"); cp ../full_initial_structure/model_${j}.pdb clean.pdb;k=$(echo "1-0.02*${i}" | bc | awk '{printf "%.2f\n", $0}'); python protein_dna_ligand_run_qdiff.py ${k};mv output.pdb output_${i}.pdb;mv output.dcd output_${i}.dcd;cd ..;done

#for i in $(seq 1 50)
#do
#  cp -r template umbrella_${i}
#  cd umbrella_${i}
#  j=$(printf "%03d" "$i")
#  cp ../full_initial_structure/model_${j}.pdb start.pdb
#  python protein_dna_ligand_build.py
#  python protein_dna_run_qinterface.py start.pdb ${j[$index]}
#  k=$(echo "0.02*${i}" | bc | awk '{printf "%.2f\n", $0}')
#  python protein_dna_ligand_run_qinterface.py ${k} >> output_energy.txt
#  cd ..
#done

