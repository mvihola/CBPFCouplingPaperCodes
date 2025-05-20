push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))

using Random, LabelledArrays, LogExpFunctions

import SimpleSMC: SMCModel
import LinearAlgebra: norm

import Base.(==)
mutable struct Float64Particle
    x::Float64
    Float64Particle() = new()
end
# Equality needed for particles to make couplings work!
(==)(x::Float64Particle, y::Float64Particle) = (x.x == y.x)

# Log-density of observation y given current log volatility lv
@inline function _svLogLikelihood(lv, y, par)
    par.normConst - .5*lv - .5*y^2*exp(-lv)
end

function svLogG(p, x, scr)
    # (12) of https://doi.org/10.1016/j.jeconom.2006.07.008
    _svLogLikelihood(x.x, scr.ys[p], scr.par)
end

@inline function _mu_x(x, p, scr)
    par = scr.par
    @inbounds y_ = scr.ys[p]
    y_cond = par.rho*par.sigma*exp(-x.x/2)*y_
    par.mu + par.phi*(x.x - par.mu) + y_cond
end

# derivatives of the above wrt parameters:
diffMu_mu_x(x, p, scr) = 1 - scr.par.phi
function diffSigma_mu_x(x, p, scr) 
    @inbounds y_ = scr.ys[p]
    scr.par.rho*exp(-x.x/2)*y_
end
diffPhi_mu_x(x, p, scr) = x.x - scr.par.mu
function diffRho_mu_x(x, p, scr) 
    @inbounds y_ = scr.ys[p]
    scr.par.sigma*exp(-x.x/2)*y_
end

function svM!(x, rng, p, x_prev, scr)
    par = scr.par
    if p == 1
        x.x = randn(rng)*sqrt(par.statVar) + par.mu
    else
        # simulate from (13) of https://doi.org/10.1016/j.jeconom.2006.07.008
        y_ = scr.ys[p-1]
        sigma_p = par.sigma*sqrt(1.0-par.rho^2)
        x.x = _mu_x(x_prev, p-1, scr) + randn(rng)*sigma_p
    end
end
function svLogM(p, x_prev, x, scr)
    par = scr.par
    if p == 1
        par.normConstStat -.5*(x.x - par.mu)^2/par.statVar
    else
        # evaluate log (13) of https://doi.org/10.1016/j.jeconom.2006.07.008
        y_ = scr.ys[p-1]
        csigmarho = 1.0/(par.sigma^2*(1.0-par.rho^2))
        par.normConst2 - .5*(x.x - _mu_x(x_prev, p-1, scr))^2*csigmarho
    end
end

newReparameterised() = LVector(logitRho=0.0, logSigma=0.0, logitPhi=0.0, mu=0.0)

function parFromTheta!(par, theta, dtheta=nothing)
    # Transformed parameter values
    logisticRho = logistic(theta.logitRho)
    rho = 2logisticRho-1
    sigma = exp(theta.logSigma)
    logisticPhi = logistic(theta.logitPhi)
    phi = 2logisticPhi-1
    mu = theta.mu
    # Set also auxiliary variables:
    setParam!(par, rho, sigma, phi, mu)
    # The gradient of the transformation:
    if !isnothing(dtheta)
        dtheta.logitRho = 2logisticRho*(1-logisticRho)
        dtheta.logitPhi = 2logisticPhi*(1-logisticPhi)
        dtheta.logSigma = sigma
        dtheta.mu = 1.0
    end
    nothing
end

function thetaFromPar(par)
    theta = newReparameterised()
    theta.mu = par.mu
    theta.logSigma = log(par.sigma)
    theta.logitRho = logit((par.rho+1)/2)
    theta.logitPhi = logit((par.phi+1)/2)
    theta
end

newPar() = LVector(rho = 0.0, sigma = 0.0, phi = 0.0, mu = 0.0, statVar = 0.0, normConst = -.2log(2pi), normConst2 = 0.0, normConstStat = 0.0)

function setParam!(par, rho, sigma, phi, mu)
    par.rho = rho
    par.sigma = sigma
    par.phi = phi
    par.mu = mu
    par.normConst2 = par.normConst -.5log(sigma^2*(1-rho^2))
    par.statVar = sigma^2/(1-phi^2)
    par.normConstStat = par.normConst -.5*log(par.statVar)
    nothing
end

function setupSvModel(y_or_n::Union{Integer,Vector}; rho=-0.7, sigma=0.1, phi=0.9, mu=0.1, lookahead=false)
    par = newPar()
    setParam!(par, rho, sigma, phi, mu)
    if typeof(y_or_n) <: Integer
        n = y_or_n
        x = Float64Particle(); x_prev = Float64Particle()
        # simulate y:
        rng = Random.GLOBAL_RNG
        y = zeros(n)
        xs = zeros(n)
        scratch = (par=par, ys=y)
        for p = 1:n
            svM!(x, rng, p, x_prev, scratch)
            y[p] = randn(rng)*exp(x.x/2)
            xs[p] = x.x
            x_prev, x = x, x_prev
        end
    else
        y = y_or_n
        n = length(y)
        scratch = (par=par, ys=y)
        xs = nothing
    end
    model = SMCModel(svM!, svLogG, n, Float64Particle, typeof(scratch))

    model, svLogM, scratch, xs
end

# Derivatives of log normal pdf with respect to μ and σ^2
gradLogNormalMu(x, mu, sigma2) = (x-mu)/sigma2
gradLogNormalSigma2(x, mu, sigma2) = - .5/sigma2 + .5(x-mu)^2/sigma2^2

function addGradSvLogM!(g, p, x_prev, x, scr)
    par = scr.par
    if p == 1
        g.mu += gradLogNormalMu(x.x, par.mu, par.statVar)
        gSigma = gradLogNormalSigma2(x.x, par.mu, par.statVar)
        g.sigma +=  gSigma*(2par.sigma/(1-par.phi^2))
        g.phi += gSigma*par.sigma^2*2par.phi/(1-par.phi^2)^2
    else
        dx = x.x - _mu_x(x_prev, p-1, scr)
        csigmarho = 1/(par.sigma^2*(1.0-par.rho^2))
        diffRho_csigmarho = csigmarho * 2par.rho/(1.0-par.rho^2)
        diffSigma_csigmarho = csigmarho * (-2/par.sigma)

        g.mu += dx*diffMu_mu_x(x_prev, p-1, scr)*csigmarho
        g.sigma += dx*diffSigma_mu_x(x_prev, p-1, scr)*csigmarho - .5*dx^2*diffSigma_csigmarho - 1.0/par.sigma
        g.rho += dx*diffRho_mu_x(x_prev, p-1, scr)*csigmarho -.5*dx^2*diffRho_csigmarho + par.rho/(1-par.rho^2)
        g.phi += dx*diffPhi_mu_x(x_prev, p-1, scr)*csigmarho
    end
end

addGradSvLogG!(g, p, x, scr) = nothing

# If the norm of vector x exceeds m, (in-place) project it to B(0,m)
function capNorm!(x, m)
    n = norm(x)
    if n > m
        x .*= (m/n)
    end
    nothing
end
