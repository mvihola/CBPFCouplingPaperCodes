using Statistics
using JLD2, CairoMakie

N = 64
T = 16801
Δ = 0.01
Ts = collect(1:T)*Δ

forwardCoupling = "ParticleMaximalCoupling"
fname_base = joinpath(@__DIR__, "..", "..", "out", "calcium_$(forwardCoupling)_$(N)")

m = 1000

elapsed = zeros(m)
couplingSites = zeros(Int, T, m)

for k = 1:m
    res = load("$(fname_base)_$(k).jld2")
    elapsed[k] = res["couplingTime"]
    couplingSites[:,k] .= res["lastCoupled"] .+ 1
end

couplingTime = mapslices(maximum, couplingSites, dims=1)[:]

Qs = mapslices(x -> quantile(x, [0.05,0.25,0.5,0.75,0.95]), 
couplingSites, dims=2)

function cm_to_pt(x)
    x/2.54*72
end
cms = (12.7, 3.5)

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

ga = f[1:2, 3:6] = GridLayout()
gcd = f[1:2, 1:2] = GridLayout()
axSites = Makie.Axis(ga[1,1], ylabel=L"\tau_t", yscale = log10, xticks = [67,68,69])
set_default_axis(axSites)
axSites2 = Makie.Axis(ga[1,2], yscale = log10, xticks=[79,80], xlabel=L"t~(s)")
set_default_axis(axSites2)
hideydecorations!(axSites2; grid=false)
axSites3 = Makie.Axis(ga[1,3:4], yscale = log10)
set_default_axis(axSites3)
hideydecorations!(axSites3; grid=false)
linkyaxes!(axSites, axSites2, axSites3)

#ind = 1:T
ind = 6700:6900
ind2 = 7850:8050
ind3 = 15500:15900
function draw_quantiles(axSites, Qs, ind)
    Makie.band!(axSites, Ts[ind], Qs[ind,1], Qs[ind,5], color=:lightgray)
    Makie.band!(axSites, Ts[ind], Qs[ind,2], Qs[ind,4], color=:gray)
    Makie.lines!(axSites, Ts[ind] ,Qs[ind,3], linewidth=linewidth, color=:black)
end
draw_quantiles(axSites, Qs, ind)
draw_quantiles(axSites2, Qs, ind2)
draw_quantiles(axSites3, Qs, ind3)

axHist = Makie.Axis(gcd[1,1], xlabel=L"\tau", ylabel=L"\mathrm{Count}")

Makie.hist!(axHist, couplingTime, bins=50)
Makie.xlims!(axHist, (0,4000))
Makie.vlines!(axHist, quantile(couplingTime, [0.5,0.9,0.95,0.99]), linewidth=linewidth, color=(:red, 0.5))
colgap!(ga, 5)

save("cai3_couplings.pdf", f; pt_per_unit=1)