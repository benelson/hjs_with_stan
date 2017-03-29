#PBS -l nodes=1:ppn=10
#PBS -l walltime=02:00:00:00
#PBS -A b1011
#PBS -q astro
#PBS -N roche

cd $PBS_O_WORKDIR
module load python
python run_pystan_2comp_2xl_overlap.py
