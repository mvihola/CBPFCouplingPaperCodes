#!/bin/sh

# Run the experiment in Slurm cluster (puhti.csc.fi):
# module load julia
# sh lgEstimatesCommandlist.sh > cmds.txt
# sbatch_commandlist -t 8:00:00 -mem 4000 -commands cmds.txt

# We run each one single-thread
CMD="julia -t 1"
SCRIPT="paper-numerics/runLgEstimates.jl"

ar=0.9
mutVar=1.0
replications=1000
output_base="out/lgEstimates"

for fwd in IndexMaximalCoupling JointIndexCoupling ParticleMaximalCoupling FilterStateMaximalCoupling; do
    for T in 512 1024 2048 4096 8192 16384; do 
        for N in 2 4 8 16 32 64 128; do
          echo $CMD $SCRIPT -f $fwd -N $N \
          --ar $ar \
          -T $T \
          --replications $replications \
          --mutVar $mutVar \
          --output ${output_base}_${ar}_${mutVar}_${fwd}_${T}_${N}.jld2 
        done
    done
done
