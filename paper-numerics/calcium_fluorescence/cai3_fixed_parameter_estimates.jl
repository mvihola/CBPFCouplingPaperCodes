using Statistics

# moving_average:
include("common.jl")
include("cai3_read.jl")

# Read data
F, N, Δ, Ts = cai3_read(3, 6)

# Poisson parameter estimate:
p_est = sum(N)/Ts[end]

# Observation noise from a silent period & smoothing
N_all = findall(N .> 0)
selection = N_all[3]:N_all[4]

F_sel = F[selection]

# Moving average with window of about 100ms:
F_smth = moving_average(F_sel, 25)

# Residual looks "roughly normal":
res = F_sel - F_smth

# Estimate σ_F by standard deviation:
σ_F_est = std(res)

# Visualise:
using Plots
p1 = plot(Ts, F, label="Full data")
p2 = plot(Ts[selection], F_sel, label="Selected data")
plot!(p2, Ts[selection], F_smth, label="Smoothed")
p3 = plot(Ts[selection], res, label="Residual")
plot(p1, p2, p3, layout=(3,1))

