using CSV, DataFrames, Dates, Statistics, Plots, LaTeXStrings, Random, JLD2, RNGPool

setRNGs(12345)
Random.seed!(12345)

include("svModelGradient.jl")

using MaximalCSMC

# Load data (data 'msci' from R package 'AER')
fname = joinpath(@__DIR__,"msci_switzerland.csv")
df = CSV.read(fname, DataFrame)

model, lM, scratch, xs = setupSvModel(diff(log.(Vector(df.msci))))
(model.pScratch)() = let scratch = scratch
    scratch
end

function initializeParameters(scratch)
    par = scratch.par
    par.rho = 0.0
    par.sigma = 1.0
    par.phi = 0.0
    par.mu = log(var(scratch.ys))
    theta0 = thetaFromPar(par)
    par0 = LVector(rho=par.rho, sigma=par.sigma, phi=par.phi, mu=par.mu)
    theta0, par0
end

h(path) = path[1].x

function saveFields(s; args=(:L, :elapsedTime, :k, :ell, :maxit, :couplingTimes, :thetas, :iterations, :n))
    NamedTuple{args}(getfield(s, i) for i in args)
end