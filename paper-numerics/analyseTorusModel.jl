using JLD2, Statistics

fwdCouplings = ("ParticleMaximalCoupling", "FilterStateMaximalCoupling", "IndexMaximalCoupling", "JointIndexCoupling",)
fwdAlias = ("IMC", "JMC", "IIC", "JIC")
cost = (N -> N^2, N -> N^2, N -> N, N -> N)
Ns = [2, 4, 8, 16, 32, 64, 128]
Ts = [512, 1024, 2048, 4096, 8192, 16384]

function read_data(scenario)
    path_prefix = joinpath(@__DIR__, "..", "out", scenario)
    iterations = repeat([NaN], length(Ns), length(Ts), length(fwdCouplings))
    costs = deepcopy(iterations)
    times = deepcopy(iterations)
    for (i, N) in enumerate(Ns)
        for (j, T) in enumerate(Ts)
            for (k, fwd) in enumerate(fwdCouplings)
                fname = "$(path_prefix)_$(fwd)_$(T)_$(N).jld2"
                if isfile(fname)
                    res = load(fname)
                    iterations[i,j,k] = mean(res["results"][1])
                    costs[i,j,k] = iterations[i,j,k]*cost[k](N)
                    times[i,j,k] = res["t_elapsed"]
                end
            end
        end
    end
    iterations, costs, times
end

include("plotFunctions.jl")

I = Array{Any}(missing, 3, 3)
C = deepcopy(I); T = deepcopy(I)
as = reverse((0.1, 0.3, 0.5)); bs = reverse((0.1, 0.3, 0.5))
for i = 1:3
    for j = 1:3
        I[i,j], C[i,j], T[i,j] = read_data("torusModel_0.2_$(as[i])_$(bs[j])")
    end
end

cms = (12.7,8)
f = Figure(
    size=cm_to_pt.(cms), 
    pt_per_unit=1, 
    fontsize=8, 
    figure_padding = (1,2,1,1)
)
plot_configurations!(f, I, 3, as, bs; title=fwdAlias[3], linear_helper=true)
plot_configurations!(f, I, 4, as, bs; offset=3, hideyaxis=true, title=fwdAlias[4], legend=7, linear_helper=true)
colgap!(f.layout, 5); rowgap!(f.layout, 5)
save("torusModel_IIC_JIC.pdf", f; pt_per_unit=1)

f = Figure(
    size=cm_to_pt.(cms), 
    pt_per_unit=1, 
    fontsize=8, 
    figure_padding = (1,2,1,1)
)
plot_configurations!(f, I, 1, as, bs; title=fwdAlias[1])
plot_configurations!(f, I, 2, as, bs; offset=3, hideyaxis=true, title=fwdAlias[2], legend=true)
colgap!(f.layout, 5); rowgap!(f.layout, 5)
save("torusModel_IMC_JMC.pdf", f; pt_per_unit=1)


f = Figure(
    size=cm_to_pt.(cms), 
    pt_per_unit=1, 
    fontsize=8, 
    figure_padding = (1,2,1,1)
)
plot_configurations!(f, C, 3, as, bs; title=fwdAlias[3], ylabel="cost")
plot_configurations!(f, C, 4, as, bs; offset=3, hideyaxis=true, title=fwdAlias[4], legend=7, ylabel="cost")
colgap!(f.layout, 5); rowgap!(f.layout, 5)
save("torusModel_IIC_JIC_cost.pdf", f; pt_per_unit=1)

f = Figure(
    size=cm_to_pt.(cms), 
    pt_per_unit=1, 
    fontsize=8, 
    figure_padding = (1,2,1,1)
)
plot_configurations!(f, C, 1, as, bs; title=fwdAlias[1], ylabel="cost")
plot_configurations!(f, C, 2, as, bs; offset=3, hideyaxis=true, title=fwdAlias[2], legend=7, ylabel="cost")
colgap!(f.layout, 5); rowgap!(f.layout, 5)
save("torusModel_IMC_JMC_cost.pdf", f; pt_per_unit=1)


cms = (11.5,7)
f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=8, 
    figure_padding = (1,1,2,1)
)

plot_iterations_makie!(f, 1, I[3,1]; legend=false, ylabel=latexstring("iter. \$b=$(bs[1])\$"), showxaxes=false)
plot_iterations_makie!(f, 2, I[3,2]; legend=false, ylabel=latexstring("iter. \$b=$(bs[2])\$"), showtitles=false, showxaxes=false)
plot_iterations_makie!(f, 3, I[3,3]; legend=true, ylabel=latexstring("iter. \$b=$(bs[3])\$"), showxaxes=true, showtitles=false)
colgap!(f.layout, 5); rowgap!(f.layout, 5)
save("torusIterations_a0.1.pdf", f; pt_per_unit=1)

f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=8, 
    figure_padding = (1,1,2,1)
)
plot_iterations_makie!(f, 1, C[3,1]; legend=false, ylabel=latexstring("cost \$b=$(bs[1])\$"), showxaxes=false)
plot_iterations_makie!(f, 2, C[3,2]; legend=false, ylabel=latexstring("cost \$b=$(bs[2])\$"), showtitles=false, showxaxes=false)
plot_iterations_makie!(f, 3, C[3,3]; legend=true, ylabel=latexstring("cost \$b=$(bs[3])\$"), showxaxes=true, showtitles=false)
colgap!(f.layout, 5); rowgap!(f.layout, 5)

axs = [content(f[i,j]) for i in 1:3, j in 1:4]
linkyaxes!(axs...)

save("torusCost_a0.1.pdf", f; pt_per_unit=1)
