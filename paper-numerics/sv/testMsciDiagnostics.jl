using JLD2, LogExpFunctions, Statistics, Random, RNGPool, LaTeXStrings, StatsBase

include("testMsci.jl")

s_IMC = jldopen("MsciIMC_ADAM.jld2")["s"]
k = s_IMC.k # Use this for burn-in & thinning
N = 16

function average_theta_last_q(thetas; q=0.9)
    iterations = length(thetas)
    start = Int(floor(q*iterations))
    ind = start:iterations
    rho = mean([2logistic(thetas[i].logitRho)-1 for i=ind])
    phi = mean([2logistic(thetas[i].logitPhi)-1 for i=ind])
    mu = mean([thetas[i].mu for i=ind])
    sigma = mean([exp(thetas[i].logSigma) for i=ind])
    rho, phi, mu, sigma
end

rho, phi, mu, sigma = average_theta_last_q(s_IMC.thetas)

setParam!(scratch.par, rho, sigma, phi, mu)

Random.seed!(12345)
setRNGs(12345)

using SimpleSMC

N = s_IMC.k
io =  SMCIO{model.particle, model.pScratch}(N, model.maxn)
smc!(model, io)
ref = [model.particle() for i = 1:model.maxn]
SimpleSMC.pickParticleBS!(ref, io, lM)

# Volatility estimates based on this many 'nearly independent' samples
m = 1000

logVol = zeros(m, model.maxn)
for i = 1:m
    for j = 1:k
        csmc!(model, lM, io, ref, ref)
    end
    for j = 1:model.maxn
        logVol[i,j] = ref[j].x
    end
end
# Quantiles for plotting log-volatilities
qs = mapslices(x -> quantile(x, [0.025,0.5,0.975]), logVol, dims=1)

volSample = logVol[end,:]

using CairoMakie, Measures

function cm_to_pt(x)
    x/2.54*72
end

cms = (12.7, 4)
f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=7, 
    figure_padding = (1,1,1,1),
)

linewidth = 0.75
function set_default_axis(a; axiswidth=linewidth)
    a.spinewidth = axiswidth
    a.xtickwidth = axiswidth
    a.ytickwidth = axiswidth
    a.xgridwidth = axiswidth
    a.ygridwidth = axiswidth
end

ga = f[1, 1:2] = GridLayout()
gb = f[2, 1:2] = GridLayout()
gcd = f[1:2, 3] = GridLayout()
axMsci = Axis(ga[1,1], ylabel=L"\mathrm{MSCI}", title=L"\mathrm{Data\ and\ estimated\ volatility}")
set_default_axis(axMsci)
Makie.lines!(axMsci, df.date, df.msci, linewidth=linewidth)
hidexdecorations!(axMsci; grid=false)

yrs = collect(Date("1995-01-01"):Year(2):Date("2013-01-01"))

axVol = Axis(gb[1,1], ylabel=L"\mathrm{Log\ volatility}", yticks=collect(-12:2:-6), xlabel=L"\mathrm{Year}", xticks=Dates.value.(yrs), xtickformat=vals -> @. string(Int(year(vals))))
set_default_axis(axVol)
Makie.band!(axVol, Dates.value.(df.date[2:end]), qs[1,:], qs[3,:], color=:gray)
Makie.lines!(axVol, Dates.value.(df.date[2:end]), qs[2,:], linewidth=0.1, color=:black)
linkxaxes!(axMsci, axVol)

rowgap!(f.layout, 5)



function calculate_epsilons(y, x)
    y .* exp.(-x./2)
end

e = calculate_epsilons(diff(log.(df.msci)), volSample)
Makie.lines(e) # looks like white noise
function acf(x; beta=1.96, ax=nothing)
    if isnothing(ax)
        f = Figure()
        ax = Axis(f[1,1])
    end
    n = length(x)
    Makie.stem!(ax, autocor(x))
    Makie.hlines!(ax, [-beta/sqrt(n), beta/sqrt(n)], color=:red)
    f
end
acf(e) # Looks like white noise

# epsilon Q-Q:
axQq = Axis(gcd[1,1], title=L"\mathrm{Normal~quantiles~of~}(\epsilon_t)")
set_default_axis(axQq)
Makie.ablines!(axQq, 0, 1, linewidth=linewidth, color=:black)
Makie.qqnorm!(axQq, e, linestyle=:solid, markersize=2)

save("MsciDiagnostics.pdf", f; pt_per_unit=1)