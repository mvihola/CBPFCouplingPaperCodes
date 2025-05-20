# Joint index coupling (JIC)

struct JointIndexCoupling
end

# log of multinomial probability mass: ∑_i log(w[a_i]) 
function _logUnnormalisedMultinomial(a, w, start=0)
    L = 0.0
    for i = start+1:length(a)
        L += log(w[a[i]])
    end
    L
end

function _jointCoupledResample!(ccsmcio)
    rng = getRNG()

    N = ccsmcio.N
    ws1 = ccsmcio.smcio1.ws
    as1 = ccsmcio.smcio1.internal.as
    ws2 = ccsmcio.smcio2.ws
    as2 = ccsmcio.smcio2.internal.as
    scratch1 = ccsmcio.smcio1.internal.scratch1
    scratch2 = ccsmcio.smcio1.internal.scratch2
  
    _sampleCategoricalSorted!(as1, ws1, N-1, scratch1, scratch2, 1, rng)
    logp = _logUnnormalisedMultinomial(as1, ws1, 1)
    logq = _logUnnormalisedMultinomial(as1, ws2, 1)
    if rand(rng) <= exp(logq - logp)
        as2 .= as1
        return
    else
        while true
            _sampleCategoricalSorted!(as2, ws2, N-1, scratch1, scratch2, 1, rng)
            logp = _logUnnormalisedMultinomial(as2, ws1, 1)
            logq = _logUnnormalisedMultinomial(as2, ws2, 1)
            if rand(rng) > exp(logp - logq)
                return
            end
        end
    end
end

function _forwardCoupling!(::JointIndexCoupling, ccsmcio, M!, p, lM)

    zetas1 = ccsmcio.smcio1.zetas
    zetaAncs1 = ccsmcio.smcio1.internal.zetaAncs
    ws1 = ccsmcio.smcio1.ws
    as1 = ccsmcio.smcio1.internal.as
    pScratch = ccsmcio.smcio1.internal.particleScratch
    xref1 = ccsmcio.ref1[p]
  
    zetas2 = ccsmcio.smcio2.zetas
    zetaAncs2 = ccsmcio.smcio2.internal.zetaAncs
    ws2 = ccsmcio.smcio2.ws
    as2 = ccsmcio.smcio2.internal.as
    xref2 = ccsmcio.ref2[p]
  
    if p > 1
        _jointCoupledResample!(ccsmcio)
        _copyParticles!(zetaAncs1, zetas1, as1)
        _copyParticles!(zetaAncs2, zetas2, as2)
    end 
    @inbounds particleCopy!(zetas1[1], xref1)
    @inbounds particleCopy!(zetas2[1], xref2)
    _indexCoupledMutateParticles!(zetas1, zetas2, M!, p, zetaAncs1,
    zetaAncs2, pScratch, xref1, xref2)

end

