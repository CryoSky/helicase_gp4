#!/bin/bash
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --job-name=check
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=sj52@rice.edu
#SBATCH --mail-type=ALL
#SBATCH -o trial.out
#SBATCH -e trial.err


date
echo My job ran on: $SLURM_NODELIST
echo My jobid was $SLURM_JOBID
echo My submission directory was $SLURM_SUBMIT_DIR

echo `which python`

cd run_job
for i in $(seq 1 50);do sed -i "/OPENMM/d" output_${i}.pdb;done
for i in $(seq 1 50);do python /home/sj52/script/helicase/calculate_qdiff_in_math.py output_${i}.pdb ../model_001_protein_cg_correct.pdb ../model_050_protein_cg_correct.pdb 0.827 qdiff_af_${i}_scaled.txt;done
for i in $(seq 1 50);do python /home/sj52/script/helicase/calculate_com_two_chain.py output_${i}.pdb A F com_af_${i}.txt;done
for i in $(seq 1 50);do python /home/sj52/script/helicase/calculate_com_two_chain.py output_${i}.pdb E F com_ef_${i}.txt;done
cd ..

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

