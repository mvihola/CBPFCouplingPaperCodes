# CBPFCouplingPaperCodes

This repository provides an implementation of coupled conditional backward sampling particle filter (CCBPF) algorithms as described in: 

> Joona Karjalainen, Anthony Lee, Sumeetpal S. Singh and Matti Vihola [Mixing time of the conditional backward sampling particle filter](https://arxiv.org/abs/2312.17572), arXiv:2312.17572, 2023.


The codes run (at least) on Julia v1.11. Before running the code, you need to install the dependencies, for instance by running the following:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/awllee/SimpleSMC.jl.git"))
for p in ("RNGPool", "NonUniformRandomVariateGeneration",
          "SMCExamples", "SequentialMonteCarlo",
          "LogExpFunctions", "JLD2", "StatsBase", "CSV",
          "ProgressMeter", "DataFrames", "ArgParse",
          "LaTeXStrings", "Measures", "CairoMakie", 
          "Distributions", "Enzyme", "ComponentArrays", "MAT")
    Pkg.add(p)
end
```

The source codes are in the following two folders:

* [src](src): The generic implementation of coupled conditional backward sampling particle filter and the related stochastic gradient maximum likelihood estimation
* [paper-numerics](paper-numerics): The codes (and scripts) for reproducing the experiments in the paper; see [paper-numerics/README.md](paper-numerics/README.md).
