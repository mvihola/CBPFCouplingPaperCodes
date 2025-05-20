using RNGPool, Statistics, StatsBase, CairoMakie, JLD2, Random, LaTeXStrings

setRNGs(1234567)
Random.seed!(1234567)

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("demo.jl")
include("lgModel.jl")

ar=0.99
mutVar=0.105
obsVar=10.0
T = 2^10

# This is the most challenging scenario of lgEstimates3
model, lM, ko = setupLG0Model(T, ar, obsVar, mutVar)

runSimulations = true
if runSimulations
    
    # This is "stable":
    N1 = 64
    MF1, PF1 = visualizeCCSMCgen(model, lM, N1; forwardCoupling=ParticleMaximalCoupling(), maxit=8192)

    # This is "critical":
    N2 = 2
    MF2, PF2 = visualizeCCSMCgen(model, lM, N2; forwardCoupling=ParticleMaximalCoupling(), maxit=8192)

    # This is "stable":
    N3 = 256
    MF3, PF3 = visualizeCCSMCgen(model, lM, N3; forwardCoupling=FilterStateMaximalCoupling(N3), maxit=8192)

    # This "never critical"
    N4 = 8
    MF4, PF4 = visualizeCCSMCgen(model, lM, N4; forwardCoupling=FilterStateMaximalCoupling(N4), maxit=8192)

    # This "never critical"
    N5 = 16
    MF5, PF5 = visualizeCCSMCgen(model, lM, N5; forwardCoupling=IndexMaximalCoupling(), maxit=8192)
    
    @save "lgDemoTraces.jld2" MF1 MF2 MF3 MF4 MF5
else
    (MF1, MF2, MF3, MF4, MF5) = load("lgDemoTraces.jld2", "MF1", "MF2", "MF3", "MF4", "MF5")
end

function cm_to_pt(x)
    x/2.54*72
end

cms = (12.7,3.6) # Figure size (in cm)
f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=8, 
    figure_padding = (1,4,1,1),
)

axiswidth = 0.5

cmap = Reverse(:binary)
ax1, hm1 = Makie.heatmap(f[1,1], MF1, colormap=cmap)

ax1.title = LaTeXString("IMC \$N=$(N1-1)\$")
ax1.ylabel = L"time index $t$"
ax1.xlabel = L"iteration $k$"
ax1.spinewidth = axiswidth
ax1.xtickwidth = axiswidth
ax1.ytickwidth = axiswidth
ax1.xgridwidth = axiswidth
ax1.ygridwidth = axiswidth
ax1.xticks = [2,4,6]

ax2, hm2 = Makie.heatmap(f[1,3], MF2, colormap=cmap)
ax2.spinewidth = axiswidth
ax2.xtickwidth = axiswidth
ax2.ytickwidth = axiswidth
ax2.xgridwidth = axiswidth
ax2.ygridwidth = axiswidth
ax2.title = LaTeXString("IMC \$N=$(N2-1)\$")
hideydecorations!(ax2, grid = false)
ax2.xlabel = L"iteration $k$"
ax2.xticks = [100,200,300]

ax3, hm3 = Makie.heatmap(f[1,2], MF3, colormap=cmap)
ax3.spinewidth = axiswidth
ax3.xtickwidth = axiswidth
ax3.ytickwidth = axiswidth
ax3.xgridwidth = axiswidth
ax3.ygridwidth = axiswidth
ax3.title = LaTeXString("JMC \$N=$(N3-1)\$")
hideydecorations!(ax3, grid = false)
ax3.xlabel = L"iteration $k$"
ax3.xticks = [2,4,6,8,10]

ax4, hm4 = Makie.heatmap(f[1,4], MF4, colormap=cmap)
ax4.spinewidth = axiswidth
ax4.xtickwidth = axiswidth
ax4.ytickwidth = axiswidth
ax4.xgridwidth = axiswidth
ax4.ygridwidth = axiswidth
ax4.title = LaTeXString("JMC \$N=$(N4-1)\$")
hideydecorations!(ax4, grid = false)
ax4.xlabel = L"iteration $k$"
ax4.xticks = [100,200,300]

ax5, hm5 = Makie.heatmap(f[1,5], MF5, colormap=cmap)
ax5.spinewidth = axiswidth
ax5.xtickwidth = axiswidth
ax5.ytickwidth = axiswidth
ax5.xgridwidth = axiswidth
ax5.ygridwidth = axiswidth
ax5.title = LaTeXString("IIC \$N=$(N5-1)\$")
hideydecorations!(ax5, grid = false)
ax5.xlabel = L"iteration $k$"
ax5.xticks = [200,400]

colgap!(f.layout, 5)
save("couplingPlots.png", f, px_per_unit=8)
save("couplingPlots.pdf", f; pt_per_unit=1)



