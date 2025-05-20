using ArgParse, JLD2
using RNGPool

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
        "--width", "-w"
            help = "Width of random walk transition"
            arg_type = Float64
            default = 0.1
        "--dynPar", "-a"
            help = "Dynamic model parameter"
            arg_type = Float64
            default = 0.5
        "--obsPar", "-b"
            help = "Observation model parameter"
            arg_type = Float64
            default = 0.5
        "--maxit", "-i"
            help = "Maximum number of iterations"
            arg_type = Int
            default = typemax(Int)
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
maxit = parsed_args["maxit"]

w = parsed_args["width"]
a = parsed_args["dynPar"]
b = parsed_args["obsPar"]

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

include("torusModel.jl")

model, lM = setupTorusModel(T, w, a, b)

h(path) = Float64(path[end].x)

t_elapsed = @elapsed results = MaximalCSMC.unbiasedEstimates(model, lM, h, N, T, forwardCoupling, 1, m, true, maxit, false)

jldsave(out; results, t_elapsed)
