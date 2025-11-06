function couplingSites(model, lM, ccsmcio; independentInitialization = true,
  maxit = 1000, forwardCoupling = ParticleMaximalCoupling()) 
  
  initializeCCSMC(model, lM, ccsmcio, independentInitialization)

  couplingIndicators = BitArray(undef, model.maxn, maxit)
  couplingIndicators .= 1
  couplingInd = BitArray(undef, model.maxn)

  coupled = _couplingSites!(couplingIndicators, couplingInd, model, lM, ccsmcio, maxit, forwardCoupling)

  function findlastOrZero(A)
    result = findlast(A)
    isnothing(result) ? 0 : result
  end


  lastCoupled = mapslices(findlastOrZero, .!couplingIndicators, dims=2)[:]

  lastCoupled, coupled
end

function _couplingSites!(couplingIndicators, couplingInd, model, lM, ccsmcio, maxit, forwardCoupling)
  ref1 = ccsmcio.ref1
  ref2 = ccsmcio.ref2
  t_ = time()
  for i in 1:maxit
    ccMCBSpf!(model, ccsmcio, lM; forwardCoupling = forwardCoupling)
    @. couplingInd = ref1 == ref2
    all(couplingInd) && break
    couplingIndicators[:,i] .= couplingInd
  end
  return time() - t_
end

