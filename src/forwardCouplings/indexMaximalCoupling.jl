# Independent index coupling (IIC)

struct IndexMaximalCoupling
end

function _forwardCoupling!(::IndexMaximalCoupling, ccsmcio, M!, p, lM)

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
  
    _coupledResample!(ccsmcio)
    if p > 1
        _copyParticles!(zetaAncs1, zetas1, as1)
        _copyParticles!(zetaAncs2, zetas2, as2)
    end 
#    @inbounds particleCopy!(zetas1[1], xref1)
#    @inbounds particleCopy!(zetas2[1], xref2)
    _indexCoupledMutateParticles!(zetas1, zetas2, M!, p, zetaAncs1,
    zetaAncs2, pScratch, xref1, xref2)


end

