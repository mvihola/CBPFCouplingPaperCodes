include("testMsci.jl")

Random.seed!(12345)

theta0, par0 = initializeParameters(scratch)

m = 100
N_JMC = 16
maxit = 2000
h(path) = path[1].x
forwardCoupling = FilterStateMaximalCoupling(N_JMC)
@time res = MaximalCSMC.unbiasedEstimates(model, lM, h, N_JMC, model.maxn, forwardCoupling, 1, m, false, maxit)
q_JMC = Int(round(quantile(res[1], 0.9)))
maxit_JMC = 100*q_JMC

iterations = 1200
s_JMC = SGML(theta0, par0, model, lM, parFromTheta!, addGradSvLogG!, addGradSvLogM!, scratch, N_JMC, iterations; 
L=q_JMC, maxit=maxit_JMC, showProgress=true, 
forwardCoupling = forwardCoupling, 
optimizer = MaximalCSMC.ADAM(theta0; α=0.01))
jldsave("MsciJMC_ADAM.jld2"; s=saveFields(s_JMC))

