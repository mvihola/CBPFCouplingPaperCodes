# Run LG unbiased estimates providing arguments in command line

using ArgParse, JLD2
using RNGPool
import Statistics: mean, var

push!(LOAD_PATH, joinpath("..", @__DIR__, "..", "src"))
using MaximalCSMC

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
            default = 16
        "--replications", "-m"
            help = "Number of replications"
            arg_type = Int
            default = 100
        "--timeHorizon", "-T"
            help = "Time horizon"
            arg_type = Int
            default = 100
        "--ar", "-p"
            help = "AR(1) coefficient"
            arg_type = Float64
            default = 0.9
        "--obsVar", "-v"
            help = "Observation noise variance"
            arg_type = Float64
            default = 1.0
        "--mutVar", "-s"
            help = "Mutation noise variance"
            arg_type = Float64
            default = 1.0
        "--maxit", "-i"
            help = "Maximum number of iterations"
            arg_type = Int
            default = typemax(Int)
        "--function", "-c"
            help = "Test function: 'mean' or 'square'"
            arg_type = String
            default = "mean"
        "--output", "-o"
            arg_type = String
            help = "File name where output is stored"
            required = true
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
N = parsed_args["numParticles"]
T = parsed_args["timeHorizon"]
m = parsed_args["replications"]
out = parsed_args["output"]
fwd = parsed_args["forwardCoupling"]
ar = parsed_args["ar"]
obsVar = parsed_args["obsVar"]
mutVar = parsed_args["mutVar"]
maxit = parsed_args["maxit"]
test_function = parsed_args["function"]

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

# Reproducible seeding by command line args:
cmd = [PROGRAM_FILE, ARGS...]
setRNGs(hash(cmd) % Int)

include("lgModel.jl")

model, lM, ko = setupLG0Model(T, ar, obsVar, mutVar) # obs are all zero

if test_function == "mean"
    function h(path::Vector{Float64Particle})
        return path[1].x
    end
    h_mean = ko.smoothingMeans[1]
elseif test_function == "square"
    function h(path::Vector{Float64Particle})
        return path[1].x^2
    end
    h_mean = ko.smoothingVariances[1] + ko.smoothingMeans[1]^2
else
    error("No such test function: $(test_function)")
end

b = 1

t_elapsed = @elapsed results = MaximalCSMC.unbiasedEstimates(model, lM, h, N, T, forwardCoupling, b, m, true, maxit, false)

jldsave(out; results, h_mean, t_elapsed)
