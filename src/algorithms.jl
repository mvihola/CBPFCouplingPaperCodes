using SimpleSMC
using ProgressMeter

# initialize the coupled CSMC state such that smcio1 is L steps ahead of smcio2
function initializeCCSMC(model::SMCModel, lM::F, ccsmcio::CCSMCIO{Particle},
  independent::Bool = false, L::Int64=1) where
  {F<:Function, Particle}
  ref1 = ccsmcio.ref1
  ref2 = ccsmcio.ref2

  smcio1 = ccsmcio.smcio1
  smcio2 = ccsmcio.smcio2

  SimpleSMC.smc!(model, smcio1)
  SimpleSMC.pickParticle!(ref1, smcio1)

  if independent
    SimpleSMC.smc!(model, smcio2)
    SimpleSMC.pickParticle!(ref2, smcio2)
  else
    SimpleSMC._copyParticles!(ref2, ref1)
  end
  for _ in 1:L
    # Advance the other chain L steps
    SimpleSMC.csmc!(model, lM, smcio1, ref1, ref1)
  end
end

# coupled conditional "maximal coupling" backward-sampling particle filter
# accommodates different forward coupling methods by multiple dispatch
function ccMCBSpf!(model::SMCModel, ccsmcio::CCSMCIO{Particle}, lM::F; forwardCoupling=ParticleMaximalCoupling()) where
  {Particle, F<:Function}
  
  zetas1 = ccsmcio.smcio1.zetas
  zetaAncs1 = ccsmcio.smcio1.internal.zetaAncs
  lws1 = ccsmcio.smcio1.internal.lws
  prevLws1 = ccsmcio.prevLws1
  tmpLws1 = ccsmcio.rws1
  ws1 = ccsmcio.smcio1.ws
  as1 = ccsmcio.smcio1.internal.as
  pScratch = ccsmcio.smcio1.internal.particleScratch
  logZhats1 = ccsmcio.smcio1.logZhats
  lZ1 = 0.0

  zetas2 = ccsmcio.smcio2.zetas
  zetaAncs2 = ccsmcio.smcio2.internal.zetaAncs
  lws2 = ccsmcio.smcio2.internal.lws
  prevLws2 = ccsmcio.prevLws2
  tmpLws2 = ccsmcio.rws2
  ws2 = ccsmcio.smcio2.ws
  as2 = ccsmcio.smcio2.internal.as
  logZhats2 = ccsmcio.smcio2.logZhats
  lZ2 = 0.0
  
  engine = getRNG()

  for p = 1:ccsmcio.n
    if p > 1
      _copyParticles!(zetaAncs1, zetas1)
      _copyParticles!(zetaAncs2, zetas2)
    end

    ws1 ./= sum(ws1)
    ws2 ./= sum(ws2)

    if p == 1 || checkEqual(zetaAncs1, zetaAncs2,1)
      # if all particles equal, use the trivial coupling
      xref1 = ccsmcio.ref1[p]
      xref2 = ccsmcio.ref2[p]
      if p > 1 
        SimpleSMC._cresample!(ccsmcio.smcio1)
        _copyParticles!(zetaAncs1, zetas1, as1)
      end
      SimpleSMC._mutateParticles!(zetas1, engine, p, model.M!, zetaAncs1, pScratch, xref1)
      _copyParticles!(zetas2, zetas1)
      particleCopy!(zetas2[1], xref2)
    else
      # The function will be dispatched based on forwardCoupling type:
      _forwardCoupling!(forwardCoupling, ccsmcio, model.M!, p, lM)
    end
    # _forwardCoupling!(forwardCoupling, ccsmcio, model.M!, p, lM)

    _logWeightParticles!(lws1, p, model.lG, zetas1, pScratch)
    _logWeightParticles!(lws2, p, model.lG, zetas2, pScratch)

    @inbounds logZhats1[p] = _normalise_log_weights!(ws1, lws1)
    @inbounds logZhats2[p] = _normalise_log_weights!(ws2, lws2)

    _intermediateOutput!(ccsmcio.smcio1, p)
    _intermediateOutput!(ccsmcio.smcio2, p)
  end

  _pickParticlesBS!(ccsmcio, lM, model.lG)
end
