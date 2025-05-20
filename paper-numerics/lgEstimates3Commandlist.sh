#!/bin/sh

# Run the experiment in Slurm cluster (puhti.csc.fi):
# module load julia
# sh lgEstimates3Commandlist.sh > cmds.txt
# sbatch_commandlist -t 8:00:00 -mem 4000 -commands cmds.txt

# We run each one single-thread
CMD="julia -t 1"
SCRIPT="paper-numerics/runLgEstimates.jl"

ar=0.99
mutVar=0.105
obsVar=10
replications=1000
testFunction="square"
output_base="out/lgEstimates3"

for fwd in IndexMaximalCoupling JointIndexCoupling ParticleMaximalCoupling FilterStateMaximalCoupling; do
    for T in 512 1024 2048 4096 8192 16384; do 
        for N in 2 4 8 16 32 64 128; do
          echo $CMD $SCRIPT -f $fwd -N $N \
          --ar $ar \
          -T $T \
          --replications $replications \
          --mutVar $mutVar \
          --obsVar $obsVar \
          --function $testFunction \
          --output ${output_base}_${fwd}_${T}_${N}.jld2 
        done
    done
done
