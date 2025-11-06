using JLD2, CairoMakie, Measures, LaTeXStrings, Statistics

# Paramaetrs:
include("calcium_fluorescence_gradients.jl")
@load "cai3_estimates.jld2"
par = DefaultFluorescenceParam()
parFromTheta_calciumFluorescence!(par, thetas[end])

# Data:
include("cai3_read.jl")
F, N_cai3, Δ, Ts = cai3_read(3, 6)

# CBPF Output:
include("common.jl")

N = 64
n_CBPF = 10000
@load "cai3_cbpf_64_10000.jld2"

# C_t quantiles:
Cs_qs = mapslices(x -> quantile(x, [0.05,0.5,0.95]), allCs, dims=2)
Cs_med = Cs_qs[:,2]
Cs_upper = Cs_qs[:,3]
Cs_lower = Cs_qs[:,1]
    
# Mean of jumps:
sum_window = 50
Ns_mean = mapslices(mean, allNs, dims=2)[:]

# Residual:
res = F - allCs[:,end] 
e = res/par.σ_F

function cm_to_pt(x)
    x/2.54*72
end

cms = (12.7, 4.5)
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
gc = f[3, 1:2] = GridLayout()
gcd = f[1:3, 3] = GridLayout()
axCalcium = Makie.Axis(ga[1,1], ylabel=L"F_t~&~C_t", title=L"\mathrm{Data,~estimated~calcium~and~spikes}")
set_default_axis(axCalcium)
hidexdecorations!(axCalcium; grid=false)

Makie.lines!(axCalcium, Ts, F, label=L"F_t", linewidth=linewidth)
Makie.band!(axCalcium, Ts, Cs_lower, Cs_upper, color=:black)
Makie.lines!(axCalcium, Ts , Cs_med, linewidth=linewidth, color=:black, label=L"C_t~\mathrm{(90%~CI)}")
axislegend(axCalcium; position=:lt, framevisible = false, orientation = :horizontal, padding=(5,0,0,-12), valign=:top)


axSpikes = Makie.Axis(gb[1,1]; ylabel=L"N_t")
set_default_axis(axSpikes)
hidexdecorations!(axSpikes; grid=false)


Makie.lines!(axSpikes, Ts, moving_sum(N_cai3, sum_window), linewidth=linewidth, label=L"\mathrm{Dataset}")
Makie.lines!(axSpikes, Ts, moving_sum(Ns_mean, sum_window), linewidth=linewidth, color=:black, label=L"\mathrm{Inferred~(mean)}")
axislegend(axSpikes; position=:lt, framevisible = false, orientation = :horizontal, padding=(2,0,0,-12), valign=:top)

axAcceptance = Makie.Axis(gc[1,1]; ylabel=L"\mathrm{Acc.~rate}", yticks = [0.2,0.6,1.0], xlabel=L"t~ \mathrm{(seconds)}")
set_default_axis(axAcceptance)
Makie.lines!(axAcceptance, Ts, accepted/n_CBPF, linewidth=linewidth/2, color=:black)

linkxaxes!(axCalcium, axSpikes, axAcceptance)


axQq = Makie.Axis(gcd[1,1], title=L"\mathrm{Normal~quantiles~of~}(\epsilon_{F,t})")
set_default_axis(axQq)
Makie.ablines!(axQq, 0, 1, linewidth=linewidth, color=:black)
Makie.qqnorm!(axQq, e, linestyle=:solid, markersize=2)


rowgap!(f.layout, 5)

save("cai3_inference.pdf", f; pt_per_unit=1)