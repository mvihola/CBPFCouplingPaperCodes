include("testMsci.jl")

using CairoMakie, Measures

function cm_to_pt(x)
    x/2.54*72
end

s_IMC = jldopen("MsciIMC_ADAM.jld2")["s"]
s_JMC = jldopen("MsciJMC_ADAM.jld2")["s"]
s_IIC = jldopen("MsciIIC_ADAM.jld2")["s"]
s_JIC = jldopen("MsciJIC_ADAM.jld2")["s"]

function show_results_makie(s, i; ts=s.elapsedTime/60, old=nothing, linewidth=0.75, colormap = :seaborn_dark6, colorrange = (1,6), cms = (12.7, 7.5), axiswidth=linewidth, xmax=s_IMC.elapsedTime[end]/60)
    if isnothing(old)
        f = Figure(
            size=cm_to_pt.(cms), 
            fontsize=7, 
            figure_padding = (1,1,1,1),
        )
        function set_default_axis(a)
            a.spinewidth = axiswidth
            a.xtickwidth = axiswidth
            a.ytickwidth = axiswidth
            a.xgridwidth = axiswidth
            a.ygridwidth = axiswidth
            a.limits = (0,xmax,nothing,nothing)
        end
        ax = [Axis(f[1,1], ylabel = L"\rho", title=L"\mathrm{Estimates}")  Axis(f[1,2], ylabel= L"\mathrm{IMC}", title=L"\mathrm{Coupling\ times}");
              Axis(f[2,1], ylabel = L"\sigma")  Axis(f[2,2], ylabel= L"\mathrm{JMC}"); 
              Axis(f[3,1], ylabel = L"\phi")  Axis(f[3,2], ylabel= L"\mathrm{IIC}");
              Axis(f[4,1], ylabel = L"\mu", xlabel=L"\mathrm{Elapsed\ time\ (min)}") Axis(f[4,2], ylabel= L"\mathrm{JIC}", xlabel=L"\mathrm{Elapsed\ time\ (min)}")]
        map(set_default_axis, ax)
        map(a -> hidexdecorations!(a; grid=false), ax[1:3,:])
    else
        ax = old; f = nothing
    end
    Makie.lines!(ax[1,1], ts, color=i, colormap=colormap, colorrange=colorrange,
    [2logistic(s.thetas[i].logitRho)-1 for i=1:s.iterations], 
    linewidth=linewidth)
    Makie.lines!(ax[2,1], ts, color=i,colormap=colormap, colorrange=colorrange,
    [exp(s.thetas[i].logSigma) for i=1:s.iterations], 
    linewidth=linewidth)
    Makie.lines!(ax[3,1], ts, color=i,colormap=colormap, colorrange=colorrange, 
    [2logistic(s.thetas[i].logitPhi)-1 for i=1:s.iterations], 
    linewidth=linewidth)
    Makie.lines!(ax[4,1], ts, color=i,colormap=colormap, colorrange=colorrange,
    [s.thetas[i].mu for i=1:s.iterations], 
    linewidth=linewidth)
    Makie.lines!(ax[i,2], ts[2:end], s.couplingTimes[2:end], color=i, colormap=colormap, colorrange=colorrange, linewidth=linewidth)
    f, ax
end

f, ax = show_results_makie(s_IMC, 1)
show_results_makie(s_JMC, 2; old=ax)
show_results_makie(s_IIC, 3; old=ax)
show_results_makie(s_JIC, 4; old=ax)
colgap!(f.layout, 12); rowgap!(f.layout, 5)

save("MLE_time_and_couplingTimes.pdf", f; pt_per_unit=1)
