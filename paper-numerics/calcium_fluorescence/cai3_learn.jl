include("calcium_fluorescence_gradients.jl")
include("smcGradientDescent.jl")
include("common.jl")
include("cai3_read.jl")
using Random, JLD2, RNGPool

# Read data
F, N_cai3, Δ, Ts = cai3_read(3, 6)
T = length(F)

scr = FluorescenceScratch(F)
theta0 = thetaFromPar(scr.par)
scr.c = merge(scr.c, (Δ = Δ,))

# Set fixed parameters (see "cai3_fixed_parameter_estimates.jl")
theta0.log_p = log(0.3809021447991138)
theta0.log_σ_F = log(0.027858304707911667)

# Set up FK-model:
model = SMCModel(M_calcium!, lG_fluorescence, T, CalciumParticle, typeof(scr))

# Define dummy constructor for scratch type
(model.pScratch)() = let scratch = scr
    scratch
end

# Number of particles & SGD iterations
N = 64
iterations = 10000

Random.seed!(12345)
setRNGs(12345)

# Keep these parameters fixed:
omit_parameters = (:p, :σ_F)

s = smcSGML(theta0, model, lM_calcium, parFromTheta_calciumFluorescence!, addGradLogG_fluorescence!, addGradLogM_calcium!, N, iterations; gradientClipping! = g -> disable_gradients!(g, omit_parameters), conditional = true)

# Save the evolution of the estimates
jldsave("cai3_estimates.jld2"; thetas = s.thetas)

# Visualisation of the progress of estimation
using Plots
par_ = similar(scr.par)
Ps = hcat([(parFromTheta_calciumFluorescence!(par_, th); copy(par_)) for th in s.thetas]...)
p = plot(Ps', label=hcat(ComponentArrays.labels(par_)...))

