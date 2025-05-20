include("testMsci.jl")

Random.seed!(12345)

theta0, par0 = initializeParameters(scratch)

m = 100
N_IIC = 16
maxit = 2000
h(path) = path[1].x
forwardCoupling = IndexMaximalCoupling()
@time res = MaximalCSMC.unbiasedEstimates(model, lM, h, N_IIC, model.maxn, forwardCoupling, 1, m, false, maxit)
# Set parameters based on quantile of meeting times
q_IIC = Int(round(quantile(res[1], 0.9)))
maxit_IIC = 100*q_IIC

iterations = 75
s_IIC = SGML(theta0, par0, model, lM, parFromTheta!, addGradSvLogG!, addGradSvLogM!, scratch, N_IIC, iterations; 
L=q_IIC, maxit=maxit_IIC, showProgress=true, 
forwardCoupling = forwardCoupling, 
optimizer = MaximalCSMC.ADAM(theta0; α=0.01))
jldsave("MsciIIC_ADAM.jld2"; s=saveFields(s_IIC))

