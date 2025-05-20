include("testMsci.jl")

Random.seed!(12345)

theta0, par0 = initializeParameters(scratch)

m = 50
N_JIC = 16
maxit = 4000
h(path) = path[1].x
forwardCoupling = JointIndexCoupling()
@time res = MaximalCSMC.unbiasedEstimates(model, lM, h, N_JIC, model.maxn, forwardCoupling, 1, m, false, maxit)
# Set parameters based on quantile of meeting times
q_JIC = Int(round(quantile(res[1], 0.9)))
maxit_JIC = 100*q_JIC

iterations = 40
s_JIC = SGML(theta0, par0, model, lM, parFromTheta!, addGradSvLogG!, addGradSvLogM!, scratch, N_JIC, iterations; 
L=q_JIC, maxit=maxit_JIC, showProgress=true, 
forwardCoupling = forwardCoupling, 
optimizer = MaximalCSMC.ADAM(theta0; α=0.01))
jldsave("MsciJIC_ADAM.jld2"; s=saveFields(s_JIC))



