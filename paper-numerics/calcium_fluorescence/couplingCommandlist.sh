#!/bin/sh

# Run the experiment in Slurm cluster (puhti.csc.fi):
# module load julia
# sh couplingCommandlist.sh 1 1000 > cmds.txt
# sbatch_commandlist -t 8:00:00 -mem 4000 -A project_2001274 -commands cmds.txt

# We run each one single-thread
CMD="julia -t 1"
SCRIPT="paper-numerics/calcium_fluorescence/cai3_couplings.jl"
seed_base="12345"
output_base="out/calcium"

for fwd in ParticleMaximalCoupling; do
  for s in $(seq $1 $2); do 
    for N in 64; do
      s_=$((seed_base + s))
      echo $CMD $SCRIPT \
      -f $fwd -N $N -s $s_ \
      --output ${output_base}_${fwd}_${N}_${s}.jld2 
    done
  done    
done
