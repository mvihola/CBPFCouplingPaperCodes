include("testMsci.jl")

Random.seed!(12345)

theta0, par0 = initializeParameters(scratch)

m = 100

N_IMC = 16
maxit = 1000
forwardCoupling = ParticleMaximalCoupling()
@time res = MaximalCSMC.unbiasedEstimates(model, lM, h, N_IMC, model.maxn, forwardCoupling, 1, m, false, maxit)
q_IMC = Int(round(quantile(res[1], 0.9)))
maxit_IMC = 100*q_IMC

iterations = 4000
s_IMC = SGML(theta0, par0, model, lM, parFromTheta!, addGradSvLogG!, addGradSvLogM!, scratch, N_IMC, iterations; 
  L=q_IMC, maxit=maxit_IMC, showProgress=true, 
  forwardCoupling = forwardCoupling,
  optimizer = MaximalCSMC.ADAM(theta0; α=0.01))
jldsave("MsciIMC_ADAM.jld2"; s=saveFields(s_IMC))
