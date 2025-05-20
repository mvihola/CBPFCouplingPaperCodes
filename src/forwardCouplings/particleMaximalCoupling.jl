# Independent maximal coupling (IMC), i.e. product of maximal couplings

struct ParticleMaximalCoupling
end

function _evaluatePDF(y::Particle, ws::Vector{Float64}, xs::Vector{Particle},
    lM::F, p::Int64, pScratch::ParticleScratch) where {Particle, F<:Function,
    ParticleScratch}
    v = 0.0
    N = length(ws)
    for i in 1:N
        v += ws[i]*exp(lM(p, xs[i], y, pScratch))
    end
    return v
end

function _maximalCouple!(x::Particle, y::Particle, ws1::Vector{Float64}, xs1::Vector{Particle},
    ws2::Vector{Float64}, xs2::Vector{Particle}, M!::F1, lM::F2, p::Int64,
    pScratch::ParticleScratch) where {Particle, F1<:Function, F2<:Function,
    ParticleScratch}
    _sample!(x, ws1, xs1, M!, p, pScratch)
    if p == 1
        @inbounds particleCopy!(y, x)
        return
    end
    if any(isnan, ws1) || any(isnan, ws2)
        _sample!(y, ws2, xs2, M!, p, pScratch)
        return
    end
    pdf1 = _evaluatePDF(x, ws1, xs1, lM, p, pScratch)
    pdf2 = _evaluatePDF(x, ws2, xs2, lM, p, pScratch)
    rng = getRNG()
    if rand(rng) < pdf2/pdf1
        particleCopy!(y, x)
    else
        # count = 0
        while true
            # count += 1
            # if count > 10000
            #   error("too many iterations")
            # end
            _sample!(y, ws2, xs2, M!, p, pScratch)
            pdf1 = _evaluatePDF(y, ws1, xs1, lM, p, pScratch)
            pdf2 = _evaluatePDF(y, ws2, xs2, lM, p, pScratch)
            if rand(rng) > pdf1/pdf2
                return
            end
        end
    end
end

function _maximallyCoupledMutateParticles!(zetas1::Vector{Particle},
    zetas2::Vector{Particle}, zetaAncs1::Vector{Particle},
    zetaAncs2::Vector{Particle}, ws1::Vector{Float64}, ws2::Vector{Float64},
    M!::F1, p::Int64, pScratch::ParticleScratch, lM::F2, xref1::Particle,
    xref2::Particle) where {Particle, F1<:Function, ParticleScratch, F2<:Function}
    @inbounds particleCopy!(zetas1[1], xref1)
    @inbounds particleCopy!(zetas2[1], xref2)
    rng = getRNG()
    for j in 2:length(zetas1)
        _maximalCouple!(zetas1[j], zetas2[j], ws1, zetaAncs1, ws2, zetaAncs2, M!, lM, p, pScratch)
    end
end

function _forwardCoupling!(::ParticleMaximalCoupling,  ccsmcio, M!, p, lM)

    zetas1 = ccsmcio.smcio1.zetas
    zetaAncs1 = ccsmcio.smcio1.internal.zetaAncs
    ws1 = ccsmcio.smcio1.ws
    pScratch = ccsmcio.smcio1.internal.particleScratch
    xref1 = ccsmcio.ref1[p]
  
    zetas2 = ccsmcio.smcio2.zetas
    zetaAncs2 = ccsmcio.smcio2.internal.zetaAncs
    ws2 = ccsmcio.smcio2.ws
    xref2 = ccsmcio.ref2[p]

    @inbounds particleCopy!(zetas1[1], xref1)
    @inbounds particleCopy!(zetas2[1], xref2)
    for j in 2:length(zetas1)
        _maximalCouple!(zetas1[j], zetas2[j], ws1, zetaAncs1, ws2, zetaAncs2, M!, lM, p, pScratch)
    end
end
