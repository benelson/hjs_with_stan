#PBS -l nodes=1:ppn=20
#PBS -l walltime=07:00:00:00
#PBS -A b1011
#PBS -q astro
#PBS -N roche

cd $PBS_O_WORKDIR
module load python
python run_pystan_2comp_2xl_overlap.py
