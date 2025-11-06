using ComponentArrays, LogExpFunctions

include("calcium_fluorescence_model.jl")

function parFromTheta_calciumFluorescence!(par, theta, dtheta=nothing)
    par.σ_F = exp(theta.log_σ_F)
    par.σ_C = exp(theta.log_σ_C)
    par.τ = exp(theta.log_τ)
    par.a = logistic(theta.logit_a)
    par.c_0 = theta.c_0
    par.p = exp(theta.log_p)
    if !isnothing(dtheta)
        dtheta.log_σ_F = par.σ_F
        dtheta.log_σ_C = par.σ_C
        dtheta.log_τ = par.τ
        dtheta.logit_a = 2par.a*(1.0-par.a)
        dtheta.c_0 = 1.0
        dtheta.log_p = par.p
    end
    nothing
end

function thetaFromPar(par)
    ComponentVector(
        log_σ_F = log(par.σ_F), 
        log_σ_C = log(par.σ_C),
        log_τ = log(par.τ),
        logit_a=logit(par.a), 
        c_0=0par.c_0,
        log_p=log(par.p))
end
