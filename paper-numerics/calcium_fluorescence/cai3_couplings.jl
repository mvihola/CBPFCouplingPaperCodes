using ArgParse, JLD2, RNGPool
push!(LOAD_PATH, joinpath("..", @__DIR__, "..", "..", "src"))
using MaximalCSMC
using SequentialMonteCarlo: SMCModel

include("calcium_fluorescence_reparam.jl")
include("cai3_read.jl")
include("couplingSites.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--forwardCoupling", "-f"
            help = "Forward coupling scheme ∈ {ParticleMaximalCoupling, FilterStateMaximalCoupling, IndexMaximalCoupling, JointIndexMaximalCoupling}"
            arg_type = String
            default = "ParticleMaximalCoupling"
        "--numParticles", "-N"
            help = "Number of particles"
            arg_type = Int
            default = 64
        "--maxit", "-m"
            help = "Max iterations"
            arg_type = Int
            default = 8000
        "--seed", "-s"
            help = "RNG seed"
            arg_type = Int
            required = true
        "--output", "-o"
            arg_type = String
            help = "File name where output is stored"
            required = true
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

N = parsed_args["numParticles"]
seed = parsed_args["seed"]
out = parsed_args["output"]
maxit = parsed_args["maxit"]
fwd = parsed_args["forwardCoupling"]

if fwd == "ParticleMaximalCoupling"
    forwardCoupling = ParticleMaximalCoupling()
elseif fwd == "FilterStateMaximalCoupling"
    forwardCoupling = FilterStateMaximalCoupling(N)
elseif fwd == "IndexMaximalCoupling"
    forwardCoupling = IndexMaximalCoupling()
elseif fwd == "JointIndexCoupling"
    forwardCoupling = JointIndexCoupling()
else
    error("Unknown forward coupling: $fwd")
end

# Read data
F, N_cai3, Δ, Ts = cai3_read(3, 6; prefix="$(@__DIR__)/")
T = length(F)

scr = FluorescenceScratch(F)

# Set parameters corresponding to the last estimates:
est = load("$(@__DIR__)/cai3_estimates.jld2")
theta_est = est["thetas"][end]
parFromTheta_calciumFluorescence!(scr.par, theta_est)

scr.c = merge(scr.c, (Δ = Δ, ))

# Set up FK-model:
model = SMCModel(M_calcium!, lG_fluorescence, T, CalciumParticle, typeof(scr))

# Define dummy constructor for scratch type
(model.pScratch)() = let scratch = scr
    scratch
end

ccsmcio = CCSMCIO{model.particle, model.pScratch}(N, model.maxn)

setRNGs(seed)

lastCoupled, couplingTime = couplingSites(model, lM_calcium, ccsmcio; independentInitialization = true, maxit = maxit, forwardCoupling=forwardCoupling)

jldsave(out; lastCoupled=lastCoupled, couplingTime=couplingTime)

