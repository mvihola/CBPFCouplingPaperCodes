include("calcium_fluorescence_gradients.jl")
include("common.jl")
include("cai3_read.jl")

using Random, JLD2, SequentialMonteCarlo, RNGPool

# Read data
F, N_cai3, Δ, Ts = cai3_read(3, 6)
T = length(F)

scr = FluorescenceScratch(F)
# Set parameters corresponding to the last estimates:
@load "cai3_estimates.jld2"
theta_est = thetas[end]
parFromTheta_calciumFluorescence!(scr.par, theta_est)
scr.c = merge(scr.c, (Δ = Δ, inflated_poisson=true, λ=0.03))

# Set up FK-model:
model = SMCModel(M_calcium!, lG_fluorescence, T, CalciumParticle, typeof(scr))

# Define dummy constructor for scratch type
(model.pScratch)() = let scratch = scr
    scratch
end

Random.seed!(12345)
setRNGs(12345)

# Number of particles
N = 64
# Setup IO
smcio = SMCIO{model.particle, model.pScratch}(N, T, 1, true)

# Number of CBPF iterations
n_CBPF = 10000

X_smc, accepted = run_cbpf(model, smcio, lM_calcium; n_CBPF = n_CBPF)

allCs = hcat([[x.C for x in X] for X in X_smc]...)
allNs = hcat([[x.N for x in X] for X in X_smc]...)

X_smc = nothing
GC.gc()

jldsave("cai3_cbpf_$(N)_$(n_CBPF).jld2"; allCs=allCs, allNs=allNs, accepted=accepted)
