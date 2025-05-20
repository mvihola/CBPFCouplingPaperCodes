using ProgressMeter

# Calculate ∇ log p(x_{1:n}, y_{1:n}) for one path x_{1:n} using
# ∇ log G and ∇ log M
function logModelGradient!(g, x, addGradLogG!::Function, addGradLogM!::Function, scratch)
    g .= 0.0
    @inline x_prev = x[1]
    for p in 1:length(x)
        @inline x_cur = x[p]
        addGradLogG!(g, p, x_cur, scratch)
        addGradLogM!(g, p, x_prev, x_cur, scratch)
        x_prev = x_cur
    end
    nothing
end

function logModelLikelihood(x, model, lM::Function, scratch)
  @inline x_prev = x[1]
  L = 0.0
  for p in 1:length(x)
    @inline x_cur = x[p]
    L += model.lG(p, x_cur, scratch)
    L += lM(p, x_prev, x_cur, scratch)
    x_prev = x_cur
  end
  L
end

# (3.6) of https://arxiv.org/pdf/2206.05691
v_t(t, L, k, ell) = floor((t-k)/L) - ceil(max(L, t-ell)/L) + 1

"""
   SGML(theta0, par0, model, lM, parFromTheta!, addGradLogG!, addGradLogM!, scratch, N, iterations; kwargs)

Stochastic gradient maximum likelihood estimation.

# Inputs:
  - `theta0`: Initial value for the parameter
  - `par0`: Initial value of the transformed parameters (`theta[i]` is transformed into `par[i]`)
  - `parFromTheta!`: Function `parFromTheta!(par, theta, dtheta)` sets parameters `par` corresponding to `theta` and calculates derivative of the transformation to `dtheta`
  - `addGradLogG!`: Function `addGradLogG!(g, p, x, scratch)` adds the contribution of ∇ log G_p(x) to vector `g`.
  - `addGradLogM!` Function `addGradLogG!(g, p, x_prev, x_cur, scratch)` adds the contribution of ∇ log M_p(x_{p-1}, x_p) to vector `g`.
  - `scratch`: The model scratch (including parameters in `scratch.par`)
  - `N`: Number of particles
  - `iterations`: Number of iterations

# Keyword arguments (with default values):
  - `n = model.maxn`: Time horizon
  - `L = 1`: Lag of coupled chains
  - `k = L`: `Burn-in` of coupling
  - `ell = 5k`: last averaged estimate
  - `forwardCoupling = ParticleMaximalCoupling()`: forward coupling strategy
  - `independentInitialization = false`
  - `maxit = typemax(Int64)`: maximum number of iterations allowed before coupling happens
  - `stepSize = k -> k^(-0.66)`: Step size sequence
  - `showProgress = true`
  - `gradientClipping! = (x -> nothing)`: Function that post-processes model gradients (use clipping for stability)

# Output: Named tuple including the following fields:
  - `thetas`: Values of `theta` over iterations
  - `couplingTimes`: Times before coupling (>= ell)
  - `elapsedTime`: Elapsed system time in iterations
""" 
function SGML(theta0, par0, model, lM, parFromTheta!, addGradLogG!, addGradLogM!, scratch, N, iterations; 
    n = model.maxn, # Length of time series
    L = 1,    # Lag; set as "high quantile of meeting times"
    k = L,    # heuristics from https://arxiv.org/pdf/2206.05691
    ell = 5k, # -"-
    forwardCoupling = ParticleMaximalCoupling(), 
    independentInitialization = false,
    maxit = typemax(Int64), # Max iterations for the unbiased estimate
    optimizer = ADAM(theta0), #SGD(theta0; c=1.0, γ=0.66),
    showProgress = true,
    gradientClipping! = (x -> nothing) # In-place post-processing of gradients
    )

    # Progress meter
    progress = Progress(iterations, dt=10, enabled=showProgress)

    # Initialize data structures
    ccsmcio = CCSMCIO{model.particle, model.pScratch}(N, n)

    # Gradients are calculated to these:
    dtheta = similar(theta0)
    g_tmp = similar(par0)
    g = similar(par0)
    grad_theta = Vector(par0)

    # (Transformed) parameters are stored here:
    thetas = [similar(theta0) for i=1:iterations]
    # Initialize:
    thetas[1] .= theta0 

    # Coupling & wall clock times
    couplingTimes = zeros(iterations)
    elapsedTime = zeros(iterations)

    # Dump the entire "state" into a named tuple
    s = (io=ccsmcio, thetas=thetas, dtheta=dtheta, g=g, g_tmp=g_tmp, grad_theta=grad_theta, forwardCoupling=forwardCoupling, independentInitialization=independentInitialization, maxit=maxit, optimizer=optimizer, model=model, lM=lM, parFromTheta! = parFromTheta!, addGradLogG!, addGradLogM!, scratch=scratch, n=n, L=L, k=k, ell=ell, iterations=iterations, couplingTimes=couplingTimes, elapsedTime=elapsedTime, gradientClipping! = gradientClipping!)

    _SGML_worker!(s, progress)

    s
end

function _SGML_worker!(s, progress)
  t_ = time()
  for i = 2:s.iterations
    # Set parameters corresponding to previous theta & calculate transformation derivatives
    @inbounds theta_ = s.optimizer.θ
    s.parFromTheta!(s.scratch.par, theta_, s.dtheta)

    # Calculate the unbiased gradient estimate
    @inbounds s.couplingTimes[i] = _unbiasedGradientEstimate!(s)

    # Calculate *negative* gradient:
    s.grad_theta .= -s.g
    # Account for the transformation:
    s.grad_theta .*= s.dtheta

    # Optimizer step
    step!(s.optimizer, s.grad_theta, i-1)
    
    @inbounds s.thetas[i] .= s.optimizer.θ
#    @inbounds s.thetas[i] .= theta_ .+ eta_ .* s.grad_theta

    # Update progress
    ProgressMeter.update!(progress, i)
    @inbounds s.elapsedTime[i] = time() - t_
  end
end

function _unbiasedGradientEstimate!(s)
  # Initialize with given lag L
  initializeCCSMC(s.model, s.lM, s.io, s.independentInitialization, s.L)

  ref1 = s.io.ref1
  ref2 = s.io.ref2

  s.g .= 0.0

  # Common weight for every gradient value:
  w = 1.0/(s.ell - s.k + 1)
  coupled = false
  couplingTime = s.maxit

  for t = 1:s.maxit
    ccMCBSpf!(s.model, s.io, s.lM; forwardCoupling=s.forwardCoupling)

    if !coupled && MaximalCSMC.checkEqual(ref1, ref2)
      coupled = true
      couplingTime = t
    end

    # The "MCMC" average part
    if s.k <= t <= s.ell
      logModelGradient!(s.g_tmp, ref1, s.addGradLogG!, s.addGradLogM!, s.scratch)
      s.gradientClipping!(s.g_tmp)
      s.g .+= w .* s.g_tmp
    end

    # The "Bias correction" part:
    if !coupled && t >= s.k + s.L
      w_ = v_t(t, s.L, s.k, s.ell)*w
      logModelGradient!(s.g_tmp, ref1, s.addGradLogG!, s.addGradLogM!, s.scratch)
      s.gradientClipping!(s.g_tmp)
      s.g .+= w_ .* s.g_tmp
      logModelGradient!(s.g_tmp, ref2, s.addGradLogG!, s.addGradLogM!, s.scratch)
      s.gradientClipping!(s.g_tmp)
      s.g .-= w_ .* s.g_tmp
    end

    # Stop if coupled
    if t >= s.ell && coupled
      return couplingTime
    end
  end
  return couplingTime
end