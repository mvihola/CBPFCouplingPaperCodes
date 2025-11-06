using ProgressMeter, SequentialMonteCarlo
push!(LOAD_PATH, "../../src")
using MaximalCSMC

function smcSGML(theta0, model, lM, parFromTheta!, addGradLogG!, addGradLogM!, N, iterations; 
    n = model.maxn, # Length of time series
    optimizer = MaximalCSMC.ADAM(theta0), #SGD(theta0; c=1.0, γ=0.66),
    showProgress = true,
    gradientClipping! = (x -> nothing), # In-place post-processing of gradients
    conditional = false) # Whether to use CSMC (with BS)

    # Progress meter
    progress = Progress(iterations, dt=10, enabled=showProgress)

    # Initialize data structures
    io = SMCIO{model.particle, model.pScratch}(N, n, 1, true)
    path = [model.particle() for _ in 1:T]

    scratch = model.pScratch()

    # Gradients are calculated to these:
    dtheta = similar(theta0)
    g = similar(scratch.par)
    grad_theta = similar(scratch.par)

    # (Transformed) parameters are stored here:
    thetas = [similar(theta0) for i=1:iterations]
    # Initialize:
    thetas[1] .= theta0 

    # Coupling & wall clock times
    couplingTimes = zeros(iterations)
    elapsedTime = zeros(iterations)

    # Dump the entire "state" into a named tuple
    s = (io=io, thetas=thetas, dtheta=dtheta, g=g, grad_theta=grad_theta, optimizer=optimizer, model=model, lM=lM, parFromTheta! = parFromTheta!, addGradLogG! = addGradLogG!, addGradLogM! = addGradLogM!, scratch=scratch, n=n, iterations=iterations, couplingTimes=couplingTimes, path=path, elapsedTime=elapsedTime, gradientClipping! = gradientClipping!)

    if conditional
        # Initialize:
        smc!(s.model, s.io)
        SequentialMonteCarlo.pickParticle!(s.path, s.io)
        MaximalCSMC._SGML_worker!(s, progress; gradientEstimate! = csmcGradientEstimate!)
    else
        MaximalCSMC._SGML_worker!(s, progress; gradientEstimate! = smcGradientEstimate!)
    end

    s
end

function smcGradientEstimate!(s)
    smc!(s.model, s.io)
    SequentialMonteCarlo.pickParticle!(s.path, s.io)
    MaximalCSMC.logModelGradient!(s.g, s.path, s.addGradLogG!, s.addGradLogM!, s.scratch)
    s.gradientClipping!(s.g)
    1.0
end

function csmcGradientEstimate!(s)
    csmc!(s.model, s.io, s.path, s.path)
    SequentialMonteCarlo.pickParticleBS!(s.path, s.io, s.lM)
    MaximalCSMC.logModelGradient!(s.g, s.path, s.addGradLogG!, s.addGradLogM!, s.scratch)
    s.gradientClipping!(s.g)
    1.0
end

function disable_gradients!(g, fnames)
    for fname in fnames
        g[fname] = 0.0
    end
    nothing
end
