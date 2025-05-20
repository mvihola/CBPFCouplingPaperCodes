#!/bin/sh

# Run the experiment in Slurm cluster (puhti.csc.fi):
# module load julia
# sh torusCommandlist.sh > cmds.txt
# sbatch_commandlist -t 4:00:00 -mem 4000 -commands cmds.txt

# We run each one single-thread
CMD="julia -t 1"
SCRIPT="paper-numerics/runTorusModel.jl"

replications=100
w=0.2
output_base="out/torusModel"

for fwd in IndexMaximalCoupling JointIndexCoupling ParticleMaximalCoupling FilterStateMaximalCoupling; do
    for T in 512 1024 2048 4096 8192 16384; do 
        for N in 2 4 8 16 32 64 128; do
        for a in 0.5 0.3 0.1; do
        for b in 0.5 0.3 0.1; do
          echo $CMD $SCRIPT \
          -f $fwd -N $N -T $T \
          --replications $replications \
          -a $a -b $b -w $w \
          --output ${output_base}_${w}_${a}_${b}_${fwd}_${T}_${N}.jld2 
        done
        done
        done
    done
done
