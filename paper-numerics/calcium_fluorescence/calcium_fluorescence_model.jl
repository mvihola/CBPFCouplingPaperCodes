using LinearAlgebra, ComponentArrays, Distributions, Random

mutable struct CalciumParticle
    C::Float64
    N::Int64
    CalciumParticle() = new()
end

# Define comparison method
import Base.:(==)
Base.:(==)(x::CalciumParticle, y::CalciumParticle) = x.C == y.C && x.N == y.N

DefaultFluorescenceParam() = ComponentVector(
    σ_F = 0.05, 
    σ_C = 0.02, 
    τ = 1.0, 
    a = 0.25, 
    c_0 = 0.0, 
    p = 0.5)

DefaultFluorescenceConst() = (
    α = 1.0, 
    β = 0.0, 
    C_lower = 0.0, 
    C_upper = 1.0, 
    F_lower = 0.0, 
    F_upper = 1.0, 
    N_max = 15,
    Δ = 0.01, 
    λ = 0.03,
    inflated_poisson = true)

mutable struct FluorescenceScratch{T1,T2}
    par::T1
    c::T2
    y::Vector{Float64}
end

FluorescenceScratch(y) = FluorescenceScratch(DefaultFluorescenceParam(), DefaultFluorescenceConst(), y)

function simulate_calcium_fluorescence!(scratch; rng=Random.GLOBAL_RNG)
    @assert !scratch.c.inflated_poisson
    n = length(scratch.y)
    X = [CalciumParticle() for i=1:n]
    x_ = X[1]
    for k = 1:n
        M_calcium!(X[k], rng, k, x_, scratch)
        scratch.y[k] = G_fluorescence(rng, k, X[k], scratch)
        x_ = X[k]
    end
    X
end


function M_calcium!(x, rng, k, x_, scratch)
    par = scratch.par
    c = scratch.c
    λ = c.inflated_poisson ? c.λ : par.p*c.Δ
    x.N = rand(rng, truncated(Poisson(λ), nothing, c.N_max))
    if k == 1
        x.C = rand(rng, Uniform(c.C_lower, c.C_upper))
    else
        f = c.Δ/par.τ
        m_C = (1.0 - f)*x_.C + f*par.c_0 + par.a*x_.N
        sd_C = sqrt(c.Δ)*par.σ_C
        x.C = rand(rng, truncated(Normal(m_C, sd_C), c.C_lower, c.C_upper))
    end
    nothing
end

function lM_calcium_worker(k, x_, x, c, par, inflated_poisson=false)
    λ = inflated_poisson ? c.λ : par.p*c.Δ
    logpdf_N = logpdf(truncated(Poisson(λ), nothing, c.N_max), x.N)
    if k == 1
        return logpdf_N
    else
        f = c.Δ/par.τ
        m_C = (1.0 - f)*x_.C + f*par.c_0 + par.a*x_.N
        sd_C = sqrt(c.Δ)*par.σ_C
        return logpdf_N + logpdf(truncated(Normal(m_C, sd_C), c.C_lower, c.C_upper), x.C)
    end
end

function lM_calcium(k, x_, x, scratch)
    par = scratch.par
    c = scratch.c
    @inline lM_calcium_worker(k, x_, x, c, par, c.inflated_poisson)
end

# Simulate from the observation model (just for simulation!)
function G_fluorescence(rng, k, x, scratch) 
    par = scratch.par
    c = scratch.c
    rand(rng, truncated(Normal(c.α*x.C + c.β, par.σ_F), c.F_lower, c.F_upper))
end

function lG_fluorescence_worker(k, x, y_, c, par)
    logpdf(truncated(Normal(c.α*x.C + c.β, par.σ_F), c.F_lower, c.F_upper), y_)
end

function lG_fluorescence(k, x, scratch)
    @inbounds y_ = scratch.y[k]
    par = scratch.par
    c = scratch.c
    if c.inflated_poisson
        log_N = (logpdf(truncated(Poisson(par.p*c.Δ), nothing, c.N_max), x.N) - 
                 logpdf(truncated(Poisson(c.λ), nothing, c.N_max), x.N))
    else
        log_N = 0.0
    end
    log_N + @inline lG_fluorescence_worker(k, x, y_, c, par)
end

