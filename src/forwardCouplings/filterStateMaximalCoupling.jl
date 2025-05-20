# Joint maximal coupling (JMC)

using LogExpFunctions

struct FilterStateMaximalCoupling
    w::Vector{Float64}
end
FilterStateMaximalCoupling(N::Int) = FilterStateMaximalCoupling(zeros(Float64,N))

function _evaluateKernelMixtureLogPDF(W, X, logM, y, w_, p, pScratch)
    for i in eachindex(W)
        @inbounds w_[i] = log(W[i]) + logM(p, X[i], y, pScratch)
    end
    logsumexp(w_)
end

function _evaluateParticleStateLogPDF(y::Vector{Particle}, ind, ws::Vector{Float64}, xs::Vector{Particle}, lM::F, p::Int64, pScratch::ParticleScratch, w_::Vector{Float64}) where {Particle, F<:Function,
    ParticleScratch}
    v = 0.0
    for i = ind
        @inbounds v += _evaluateKernelMixtureLogPDF(ws, xs, lM, y[i], w_,  p, pScratch)
    end
    return v
end

function _filterStateMaximalCouple!(x::Vector{Particle}, y::Vector{Particle}, ws1::Vector{Float64}, xs1::Vector{Particle},
    ws2::Vector{Float64}, xs2::Vector{Particle}, M!::F1, lM::F2, p::Int64,
    pScratch::ParticleScratch, w_::Vector{Float64}) where {Particle, F1<:Function, F2<:Function,
    ParticleScratch}
    
    ind = 2:length(x)
    for i = ind
        @inbounds _sample!(x[i], ws1, xs1, M!, p, pScratch)
    end
    if p == 1
        for i = ind
            @inbounds particleCopy!(y[i], x[i])
        end
        return
    end
    if any(isnan, ws1) || any(isnan, ws2)
        for i = ind
            @inbounds _sample!(y[i], ws2, xs2, M!, p, pScratch)
        end
        return
    end
    
    logpdf1 = _evaluateParticleStateLogPDF(x, ind, ws1, xs1, lM, p, pScratch, w_)
    logpdf2 = _evaluateParticleStateLogPDF(x, ind, ws2, xs2, lM, p, pScratch, w_)
    rng = getRNG()
    if rand(rng) < exp(logpdf2 - logpdf1)
        for i = ind
            @inbounds particleCopy!(y[i], x[i])
        end
    else
        # count = 0
        while true
            # count += 1
            # if count > 10000
            #   error("too many iterations")
            # end
            for i = ind
                @inbounds _sample!(y[i], ws2, xs2, M!, p, pScratch)
            end
            logpdf1 = _evaluateParticleStateLogPDF(y, ind, ws1, xs1, lM, p, pScratch, w_)
            logpdf2 = _evaluateParticleStateLogPDF(y, ind, ws2, xs2, lM, p, pScratch, w_)
            if rand(rng) > exp(logpdf1 - logpdf2)
                return
            end
        end
    end
end

function _forwardCoupling!(fsmc::FilterStateMaximalCoupling, ccsmcio, M!, p, lM)

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

    _filterStateMaximalCouple!(zetas1, zetas2, ws1, zetaAncs1, ws2, zetaAncs2, M!, lM, p,pScratch, fsmc.w)
end
