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
    mse = deepcopy(iterations)
    times = deepcopy(iterations)
    for (i, N) in enumerate(Ns)
        for (j, T) in enumerate(Ts)
            for (k, fwd) in enumerate(fwdCouplings)
                fname = "$(path_prefix)_$(fwd)_$(T)_$(N).jld2"
                if isfile(fname)
                    res = load(fname)
                    iterations[i,j,k] = mean(res["results"][1])
                    costs[i,j,k] = iterations[i,j,k]*cost[k](N)
                    mse[i,j,k] = mean((res["results"][2].-res["h_mean"]).^2)
                    times[i,j,k] = res["t_elapsed"]
                end
            end
        end
    end
    iterations, costs, mse, times
end

include("plotFunctions.jl")

iterations, costs, mse, times = read_data("lgEstimates_0.9_1.0")
iterations2, costs2, mse2, times2 = read_data("lgEstimates2_0.99_0.105")
iterations3, costs3, mse3, times3 = read_data("lgEstimates3")


cms = (11.5,7)
f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=8, 
    figure_padding = (1,1,2,1)
)

plot_iterations_makie!(f, 1, iterations; legend=false, ylabel=L"iter., $\theta_1$", showxaxes=false)
plot_iterations_makie!(f, 2, iterations2; legend=false, ylabel=L"iter., $\theta_2$", showtitles=false, showxaxes=false)
plot_iterations_makie!(f, 3, iterations3; legend=true, ylabel=L"iter., $\theta_3$", showxaxes=true, showtitles=false)
colgap!(f.layout, 5); rowgap!(f.layout, 5)

save("lgIterations.pdf", f; pt_per_unit=1)

f = Figure(
    size=cm_to_pt.(cms), 
    fontsize=8, 
    figure_padding = (1,1,2,1)
)

plot_iterations_makie!(f, 1, costs; legend=false, ylabel=L"cost, $\theta_1$", showxaxes=false)
plot_iterations_makie!(f, 2, costs2; legend=false, ylabel=L"cost, $\theta_2$", showtitles=false, showxaxes=false)
plot_iterations_makie!(f, 3, costs3; legend=true, ylabel=L"cost, $\theta_3$", showxaxes=true, showtitles=false)
colgap!(f.layout, 5); rowgap!(f.layout, 5)

axs = [content(f[i,j]) for i in 1:3, j in 1:4]
linkyaxes!(axs...)

save("lgCost.pdf", f; pt_per_unit=1)

